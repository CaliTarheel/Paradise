"""Overview maps of a history: the polygons every so many Myr, as one sheet.

    python examples/overview.py worlds/thousand.paradise.json --every 100 --out overview.png

Each panel fills continent and ocean (by age) from the model's polygons, and draws over them, as
lines: plate outlines, boundaries by type (trenches with their teeth on the overriding side),
island and continental arcs, terranes, orogenies and hotspot trails.
"""
import argparse
import io
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from paradise import plumes                    # noqa: E402
from paradise.core import render, sphere       # noqa: E402
from paradise.history import History, Project  # noqa: E402

KINDS = dict(ridge=(210, 50, 40), trench=(40, 90, 175), transform=(40, 140, 75), converging=(200, 140, 20),
             collision=(120, 70, 30))
FEATURES = dict(island_arc=(240, 110, 40), continental_arc=(180, 40, 140), terrane=(110, 60, 160),
                lip=(125, 30, 105))
OROGENS = dict(active=(190, 40, 40), former=(140, 90, 40), old=(150, 150, 135))


def _xy(points, width):
    lon, lat = sphere.to_lonlat(np.asarray(points, float))
    return list(zip((lon + 180.) / 360. * width, (90. - lat) / 180. * (width / 2)))


def _line(pen, points, width, colour, thick=1):
    xy = _xy(points, width)
    for a, b in zip(xy[:-1], xy[1:]):
        if abs(a[0] - b[0]) < width / 2:
            pen.line([a, b], fill=colour, width=thick)


def _ring(pen, ring, width, colour, thick=1):
    _line(pen, np.vstack((ring, ring[:1])), width, colour, thick)


def panel(snapshot, age, width=720):
    world = snapshot.world
    buffer = io.BytesIO()
    render.draw(world, buffer, width=width, boundaries=False)
    image = Image.open(io.BytesIO(buffer.getvalue())).convert('RGB')
    pen = ImageDraw.Draw(image)
    for p, region in world.regions.items():
        for rings in sphere.polygons(region):
            for ring in rings:
                _ring(pen, ring, width, (20, 30, 40))
    for belt in getattr(snapshot, 'orogens', []):
        for rings in sphere.polygons(belt['shape']):
            for ring in rings:
                _ring(pen, sphere.rotate(ring, world.frame(belt['plate'])), width, OROGENS[belt['kind']], 2)
    for f in getattr(snapshot, 'features', []):
        for rings in sphere.polygons(f.shape):
            for ring in rings:
                _ring(pen, sphere.rotate(ring, world.frame(f.plate)), width, FEATURES.get(f.kind, (0, 0, 0)))
    for b in world.boundaries:
        start = 0
        for k in range(1, len(b.kinds) + 1):
            if k < len(b.kinds) and b.kinds[k] == b.kinds[start]:
                continue
            _line(pen, b.points[start:k + 1], width, KINDS.get(b.kinds[start], (0, 0, 0)), 2)
            start = k
    for t in getattr(snapshot, 'trenches', []):
        points = sphere.from_lonlat(*np.array(t['points']).T)
        left = t['overriding'] == t['plates'][0]
        for a, c in zip(points[:-1:3], points[1::3]):
            mid = sphere.unit(a + c)
            side = sphere.unit(np.cross(mid, sphere.unit(c - a)))
            tip = sphere.unit(mid + (side if left else -side) * (60. / sphere.RADIUS_KM))
            (x0, y0), (x1, y1) = _xy([mid, tip], width)
            if abs(x0 - x1) < 20:
                pen.line([(x0, y0), (x1, y1)], fill=KINDS['trench'], width=2)
    for h, plate, points, first, last in plumes.trails(getattr(snapshot, 'hotspots', [])):
        if plate in world.regions:
            _line(pen, sphere.rotate(points, world.frame(plate)), width, (255, 210, 60), 2)
    pen.rectangle([0, 0, 150, 18], fill=(255, 255, 255))
    pen.text((4, 3), '%.0f Ma: %d plates' % (age, len(world.plates())), fill=(0, 0, 0), font=ImageFont.load_default())
    return image


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('project')
    parser.add_argument('--every', type=float, default=100.)
    parser.add_argument('--width', type=int, default=720)
    parser.add_argument('--columns', type=int, default=2)
    parser.add_argument('--out', default='overview.png')
    args = parser.parse_args()
    project = Project.load(args.project)
    cache = args.project + '.cache'
    history = History(project, cache if os.path.exists(cache) else None)
    panels = [panel(s, project.age_ma(s.time_myr), args.width) for s in history.snapshots
              if abs(s.time_myr / args.every - round(s.time_myr / args.every)) < 1e-6]
    rows = (len(panels) + args.columns - 1) // args.columns
    sheet = Image.new('RGB', (args.columns * args.width, rows * args.width // 2), (255, 255, 255))
    for k, image in enumerate(panels):
        sheet.paste(image, ((k % args.columns) * args.width, (k // args.columns) * (args.width // 2)))
    sheet.save(args.out)
    print('wrote %s: %d panels' % (args.out, len(panels)))


if __name__ == '__main__':
    main()
