"""The local Paradise app: a small HTTP server for the browser editor and stepper.

    python -m paradise.server --project my-world.paradise.json --open

One project is open at a time. It is saved after every change (build, commit, going back), so the
project file is always the history so far. The browser talks to it through a few JSON routes:

    GET  /api/start     the painted start of the open project
    GET  /api/state     the world at the latest step, the map and the step list
    POST /api/build     a new project from a painting (replaces the open one's history)
    POST /api/preview   the next step with the given poles, not kept
    POST /api/commit    the next step with the given poles, kept
    POST /api/back      drop the steps after a given one
    POST /api/export    write the GPlates project beside the project file
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import os
import threading
import traceback
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

from paradise import gplates
from paradise.core import render, sphere
from paradise.history import History, Project

WEB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'web')
MAP_WIDTH = 720
TYPES = {'.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8',
         '.css': 'text/css; charset=utf-8', '.png': 'image/png'}


class App:
    def __init__(self, path):
        self.path = path
        self.lock = threading.Lock()
        self.cache = path + '.cache'
        self.history = History(Project.load(path), self.cache) if os.path.exists(path) else None

    def save(self):
        self.history.project.save(self.path)
        self.history.save_cache()

    # ------------------------------------------------------------------ views
    def start(self) -> dict:
        p = self.history.project
        return dict(name=p.name, start_age_ma=p.start_age_ma, ocean_age_myr=p.ocean_age_myr,
                    width=int(p.plates.shape[1]), height=int(p.plates.shape[0]),
                    plates=base64.b64encode(p.plates.astype('<u2').tobytes()).decode('ascii'),
                    crust=base64.b64encode(p.crust.astype('<u2').tobytes()).decode('ascii'))

    def state(self) -> dict:
        history = self.history
        project = history.project
        snapshot = history.current
        world = snapshot.world
        buffer = io.BytesIO()
        plate_of = render.draw(world, buffer, width=MAP_WIDTH, boundaries=False)
        last = project.steps[-1]['poles'] if project.steps else {}
        plates = []
        for p in world.plates():
            lon, lat = sphere.to_lonlat(_centre(world.regions[p]))
            pole = last.get(str(p))
            if pole is None:
                parent = world.parents.get(p)
                pole = last.get(str(parent), [90., 0., 0.]) if parent is not None else [90., 0., 0.]
            plates.append(dict(id=p, name=world.names.get(p, 'Plate %d' % p),
                               area_mkm2=round(sphere.area_km2(world.regions[p]) / 1e6, 2),
                               centre=[float(lon), float(lat)], pole=pole))
        steps = [dict(index=0, age_ma=project.age_ma(0.), note='Painted start')]
        for k, step in enumerate(project.steps, start=1):
            steps.append(dict(index=k, age_ma=project.age_ma(history.snapshots[k].time_myr),
                              note=step.get('note', '')))
        continent = sum(c.area_km2() for c in world.crust if c.kind == 'continent')
        return dict(name=project.name, time_myr=snapshot.time_myr, age_ma=project.age_ma(snapshot.time_myr),
                    step_myr=project.step_myr, steps=steps, plates=plates,
                    continent_share=continent / sphere.SPHERE_KM2,
                    map='data:image/png;base64,' + base64.b64encode(buffer.getvalue()).decode('ascii'),
                    plate_of=base64.b64encode(plate_of.astype('<u2').tobytes()).decode('ascii'),
                    map_width=MAP_WIDTH, lines=_lines(world, snapshot.cratons), warnings=snapshot.warnings,
                    events=_events(snapshot))

    # ------------------------------------------------------------------ actions
    def build(self, body) -> dict:
        width, height = int(body['width']), int(body['height'])
        plates = _grid(body['plates'], height, width)
        crust = _grid(body['crust'], height, width)
        project = Project(body.get('name') or 'Paradise', plates, crust,
                          start_age_ma=float(body.get('start_age_ma', 2000.)),
                          ocean_age_myr=float(body.get('ocean_age_myr', 60.)))
        self.history = History(project, self.cache)
        self.save()
        result = self.state()
        result['message'] = 'Built %d plates' % len(result['plates'])
        return result

    def preview(self, body) -> dict:
        snapshot = self.history.preview(_poles(body))
        return dict(lines=_lines(snapshot.world, snapshot.cratons), warnings=snapshot.warnings, events=_events(snapshot),
                    message='Preview of %.0f Ma' % self.history.project.age_ma(snapshot.time_myr))

    def commit(self, body) -> dict:
        self.history.commit(_poles(body), body.get('note', ''))
        self.save()
        result = self.state()
        result['message'] = 'Committed: now %.0f Ma' % result['age_ma']
        return result

    def back(self, body) -> dict:
        self.history.back_to(int(body['step']))
        self.save()
        result = self.state()
        result['message'] = 'Back at %.0f Ma' % result['age_ma']
        return result

    def export(self, body) -> dict:
        stem = os.path.splitext(os.path.splitext(self.path)[0])[0]
        directory = stem + '-gplates'
        files = gplates.export(self.history, directory)
        return dict(directory=os.path.abspath(directory), files=[os.path.basename(f) for f in files],
                    message='Exported %d files' % len(files))


def _grid(text, height, width) -> np.ndarray:
    data = np.frombuffer(base64.b64decode(text), '<u2')
    if data.size != height * width:
        raise ValueError('the painting is %d cells, not %d by %d' % (data.size, height, width))
    return data.reshape(height, width).astype(np.int64)


def _poles(body) -> dict:
    poles = {}
    for plate, pole in (body.get('poles') or {}).items():
        lat, lon, rate = (float(v) for v in pole)
        if not -90. <= lat <= 90.:
            raise ValueError('plate %s: pole latitude must be between -90 and 90' % plate)
        if rate:
            poles[int(plate)] = [lat, lon, rate]
    return poles


def _centre(region) -> np.ndarray:
    import spherely
    c = spherely.centroid(region)
    return sphere.from_lonlat(spherely.get_x(c), spherely.get_y(c))


def _lonlat(points) -> list:
    lon, lat = sphere.to_lonlat(np.asarray(points, float))
    return [[round(float(a), 3), round(float(b), 3)] for a, b in zip(lon, lat)]


def _lines(world, cratons=()) -> dict:
    """Boundaries by kind, plate and craton outlines, as lon/lat polylines for the browser to draw."""
    boundaries = []
    for b in world.boundaries:
        start = 0
        for k in range(1, len(b.kinds) + 1):
            if k == len(b.kinds) or b.kinds[k] != b.kinds[start]:
                boundaries.append(dict(kind=b.kinds[start], plates=[b.left, b.right],
                                       overriding=b.overriding[start], points=_lonlat(b.points[start:k + 1])))
                start = k
    outlines = []
    for p in world.plates():
        for rings in sphere.polygons(world.regions[p]):
            for ring in rings:
                outlines.append(dict(plate=p, points=_lonlat(np.vstack((ring, ring[:1])))))
    craton_lines = []
    for plate, shape in cratons:
        frame = world.frame(plate)
        for rings in sphere.polygons(shape):
            for ring in rings:
                world_ring = sphere.rotate(ring, frame)
                craton_lines.append(dict(plate=plate, points=_lonlat(np.vstack((world_ring, world_ring[:1])))))
    return dict(boundaries=boundaries, outlines=outlines, cratons=craton_lines)


# Bookkeeping the core does as plates move, not worth a person's attention: slivers below the
# level of detail passed to a neighbour, and pieces smaller than this given to the plate they border.
SMALL_PIECE_KM2 = 1e5


def _events(snapshot) -> list:
    """What happened during a step, once per kind and plates, with the routine bookkeeping left out."""
    events, seen = [], set()
    for entry in snapshot.entries:
        for e in entry.get('events', []):
            kind = e.get('kind')
            if kind == 'crumbs_given' or (kind == 'absorbed' and e.get('area_km2', 0.) < SMALL_PIECE_KM2):
                continue
            key = (kind, tuple(e['plates']) if 'plates' in e else e.get('plate'), e.get('overriding'))
            if key in seen:
                continue
            seen.add(key)
            events.append(dict(e, time_myr=entry.get('time_myr')))
    return events


class Handler(BaseHTTPRequestHandler):
    app: App = None

    def log_message(self, format, *args):
        pass

    def _send(self, status, body, kind='application/json'):
        data = body if isinstance(body, bytes) else json.dumps(body).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', kind)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        path = self.path.split('?')[0]
        if path in ('/api/start', '/api/state'):
            with self.app.lock:
                if self.app.history is None:
                    return self._send(200, dict(empty=True))
                return self._send(200, self.app.start() if path == '/api/start' else self.app.state())
        name = 'index.html' if path == '/' else path.lstrip('/')
        full = os.path.realpath(os.path.join(WEB, name))
        if not full.startswith(os.path.realpath(WEB) + os.sep) or not os.path.isfile(full):
            return self._send(404, dict(error='not found'))
        with open(full, 'rb') as f:
            self._send(200, f.read(), TYPES.get(os.path.splitext(full)[1], 'application/octet-stream'))

    def do_POST(self):
        routes = {'/api/build': App.build, '/api/preview': App.preview, '/api/commit': App.commit,
                  '/api/back': App.back, '/api/export': App.export}
        route = routes.get(self.path.split('?')[0])
        if route is None:
            return self._send(404, dict(error='not found'))
        length = int(self.headers.get('Content-Length') or 0)
        try:
            body = json.loads(self.rfile.read(length) or b'{}')
            with self.app.lock:
                if route is not App.build and self.app.history is None:
                    return self._send(409, dict(error='No world yet: paint one and build it.'))
                return self._send(200, route(self.app, body))
        except (ValueError, KeyError) as e:
            return self._send(400, dict(error=str(e)))
        except Exception as e:  # noqa: BLE001 - shown in the browser, with the trace in the terminal
            traceback.print_exc()
            return self._send(500, dict(error='%s: %s' % (type(e).__name__, e)))


def main(argv=None):
    parser = argparse.ArgumentParser(description='Paradise: a guided tectonic history.')
    parser.add_argument('--project', default='world.paradise.json',
                        help='project file to open, or to create when you build (default: %(default)s)')
    parser.add_argument('--port', type=int, default=8790)
    parser.add_argument('--open', action='store_true', help='open the app in the browser')
    args = parser.parse_args(argv)
    Handler.app = App(args.project)
    server = ThreadingHTTPServer(('127.0.0.1', args.port), Handler)
    url = 'http://127.0.0.1:%d/' % args.port
    print('Paradise is running at %s (project %s). Press Ctrl+C to stop.' % (url, os.path.abspath(args.project)))
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
