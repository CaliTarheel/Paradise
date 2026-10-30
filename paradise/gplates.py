"""Export a history as a GPlates project: a rotation file and feature collections, one per layer.

Model time runs forward from 0 at the painted start; GPlates counts age back from the present. Age is
`start_age_ma - time_myr`, so the start is at `start_age_ma` and the latest step at the youngest age.

Each plate's total rotation at an age is its frame then: the rotation from its own frame, which is
where its crust was at the start, to where it is. Every feature is a snapshot of one step, valid
from that step's age until the next step's, with its geometry carried back into its plate's frame.
GPlates then rotates it smoothly between steps: within a step each plate turns about one pole, and
GPlates' interpolation between two rotations about one pole is that turn.

Layers follow the hand-built Artifexia project, as far as this phase of the model goes.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pygplates

from paradise import plumes
from paradise.core import sphere

ANCHOR = 0
NOTE = 'Paradise'

BOUNDARY_TYPES = {
    'ridge': ('MidOceanRidge', 'Divergent'),
    'trench': ('SubductionZone', 'Convergent'),
    'collision': ('OrogenicBelt', 'Convergent'),
    'converging': ('UnclassifiedFeature', 'Convergent'),
    'transform': ('Transform', 'Transform'),
}


def matrix_pole(matrix) -> tuple:
    """(lat, lon, angle in degrees) of a rotation matrix."""
    matrix = np.asarray(matrix, float)
    angle = math.acos(max(-1., min(1., (np.trace(matrix) - 1.) / 2.)))
    if angle < 1e-12:
        return 90., 0., 0.
    axis = np.array([matrix[2, 1] - matrix[1, 2], matrix[0, 2] - matrix[2, 0], matrix[1, 0] - matrix[0, 1]])
    if np.linalg.norm(axis) < 1e-9:
        # Near 180 degrees: the axis is the column of (R + I) with the largest norm.
        m = matrix + np.eye(3)
        axis = m[:, int(np.argmax(np.linalg.norm(m, axis=0)))]
    lon, lat = sphere.to_lonlat(sphere.unit(axis))
    return float(lat), float(lon), math.degrees(angle)


def rotation_rows(history) -> list:
    """(plate, age, lat, lon, angle, fixed plate) rows for every plate at every sub-step it exists.

    A plate made by splitting another (a rift, a back-arc) turns relative to the plate it came from,
    as in a hand-built rotation tree, when that plate is there at every one of its rows; otherwise,
    and for every other plate, relative to the anchor plate 000.
    """
    project = history.project
    frames_at = {}   # plate -> {time: frame}
    for snapshot in history.snapshots:
        for time_myr, frames in snapshot.frames:
            for plate, frame in frames.items():
                frames_at.setdefault(plate, {}).setdefault(round(time_myr, 6), frame)
    parents = history.current.world.parents
    rows = []
    for plate, times in frames_at.items():
        parent = parents.get(plate)
        fixed = parent if parent in frames_at and all(t in frames_at[parent] for t in times) else ANCHOR
        for time_myr, frame in times.items():
            relative = frame if fixed == ANCHOR else frames_at[fixed][time_myr].T @ frame
            lat, lon, angle = matrix_pole(relative)
            rows.append((plate, project.age_ma(time_myr), lat, lon, angle, fixed))
    rows.sort(key=lambda r: (r[0], r[1]))   # GPlates reads each plate's rows youngest first
    return rows


def write_rotations(history, path):
    """A PLATES4 rotation file (.rot); each plate relative to the anchor plate 000."""
    with open(path, 'w', encoding='ascii') as f:
        for plate, age, lat, lon, angle, fixed in rotation_rows(history):
            f.write('%03d %9.4f %9.4f %10.4f %10.4f  %03d ! %s plate %d\n'
                    % (plate, age, lat, lon, angle, fixed, NOTE, plate))


def _points_on_sphere(points) -> list:
    lon, lat = sphere.to_lonlat(np.asarray(points, float))
    return list(zip(lat.tolist(), lon.tolist()))


def _polygons(geography, back) -> list:
    """pygplates polygons of a geography, carried by `back` (world to plate frame)."""
    result = []
    for rings in sphere.polygons(geography):
        if len(rings[0]) < 3:
            continue
        exterior = _points_on_sphere(sphere.rotate(rings[0], back))
        interiors = [_points_on_sphere(sphere.rotate(r, back)) for r in rings[1:] if len(r) >= 3]
        try:
            result.append(pygplates.PolygonOnSphere(exterior, interiors) if interiors
                          else pygplates.PolygonOnSphere(exterior))
        except (pygplates.InvalidPointsForPolygonConstructionError, ValueError):
            continue
    return result


def _feature(kind, geometry, name, plate, begin, end, description=None):
    feature = pygplates.Feature(pygplates.FeatureType.create_gpml(kind))
    feature.set_geometry(geometry)
    feature.set_name(name)
    feature.set_reconstruction_plate_id(int(plate))
    feature.set_valid_time(begin, end)
    if description:
        feature.set_description(description)
    return feature


def _windows(history):
    """(snapshot, begin age, end age) for every step: each shown until the next one."""
    project = history.project
    snapshots = history.snapshots
    for k, snapshot in enumerate(snapshots):
        begin = project.age_ma(snapshot.time_myr)
        if k + 1 < len(snapshots):
            end = project.age_ma(snapshots[k + 1].time_myr)
        else:
            end = begin - project.step_myr
        # GPlates takes a feature as existing from begin (older) to end (younger), inclusive at the
        # older end; nudge the younger end so consecutive snapshots do not both show at one age.
        yield k, snapshot, begin, end + 1e-3


def _hotspots(history, out):
    """Each hotspot once, fixed on the anchor plate, for as long as it lasts."""
    project = history.project
    for h in getattr(history.current, 'hotspots', []):
        lon, lat = h['point']
        end = project.age_ma(h['ended_myr']) if h['ended_myr'] is not None else -1e-3
        out['Hotspots'].append(_feature('HotSpot', pygplates.PointOnSphere(lat, lon), h['name'], ANCHOR,
                                        project.age_ma(h['born_myr']), end))


def _signature(feature) -> tuple:
    """What a feature is, but for when: its type, name, plate, description and geometry."""
    rings = []
    for geometry in feature.get_all_geometries():
        if isinstance(geometry, pygplates.PolygonOnSphere):
            parts = [geometry.get_exterior_ring_points()] + [
                geometry.get_interior_ring_points(k) for k in range(geometry.get_number_of_interior_rings())]
        elif isinstance(geometry, pygplates.PointOnSphere):
            parts = [[geometry]]
        else:
            parts = [geometry.get_points()]
        for part in parts:
            rings.append(np.round(np.array([q.to_lat_lon() for q in part], float), 5).tobytes())
    return (feature.get_feature_type().to_qualified_string(), feature.get_name(), feature.get_reconstruction_plate_id(),
            feature.get_description(), tuple(rings))


def _merged(features) -> list:
    """The features with each run of identical step snapshots written once, for the run's whole time.

    A piece of crust is the same in its plate's frame from one step to the next until something cuts
    it; written once a step, the thousand-Myr export came to 1.1 GB, against Artifexia's few MB.
    """
    result, open_runs = [], {}
    for feature in features:
        begin, end = feature.get_valid_time()
        key = _signature(feature)
        run = open_runs.get(key)
        if run is not None and abs(run.get_valid_time()[1] - begin) < .01:
            run.set_valid_time(run.get_valid_time()[0], end)
            continue
        open_runs[key] = feature
        result.append(feature)
    return result


def layers(history) -> dict:
    """{layer file name: [features]} for the whole history."""
    out = {name: [] for name in ('Plates', 'Continents', 'Cratons', 'Island arcs', 'Continental arcs',
                                 'Accreted terranes', 'Rifts', 'Failed rifts', 'Active orogenies',
                                 'Former orogenies', 'Old orogenies', 'Active LIPs', 'Former LIPs',
                                 'Hotspots', 'Hotspot trails', 'Ocean crust', 'Divergent', 'Convergent',
                                 'Transform')}
    for k, snapshot, begin, end in _windows(history):
        world = snapshot.world
        for plate in world.plates():
            back = world.frame(plate).T
            for polygon in _polygons(world.regions[plate], back):
                out['Plates'].append(_feature('UnclassifiedFeature', polygon, 'plate %d' % plate,
                                              plate, begin, end))
        for c in world.crust:
            # Crust is stored in its plate's frame already.
            for polygon in _polygons(c.shape, np.eye(3)):
                if c.kind == 'continent':
                    out['Continents'].append(_feature('ContinentalCrust', polygon, 'continent',
                                                      c.plate, begin, end,
                                                      'thickness %.1f km' % c.thickness_km))
                else:
                    born = history.project.age_ma(c.born_myr) if c.born_myr is not None else None
                    name = 'ocean crust %s Ma' % ('%.0f' % born if born is not None else 'unknown')
                    out['Ocean crust'].append(_feature('OceanicCrust', polygon, name, c.plate, begin, end))
        for plate, shape in snapshot.cratons:
            for polygon in _polygons(shape, np.eye(3)):
                out['Cratons'].append(_feature('Craton', polygon, 'craton', plate, begin, end))
        for r in getattr(snapshot, 'rifts', []):
            # The line the plate was split along, carried by the plate it was drawn on.
            geometry = pygplates.PolylineOnSphere(_points_on_sphere(r['points']))
            # A rift whose two sides are one plate again (welded back) failed.
            failed = r['plate'] not in world.regions or r['new'] not in world.regions
            out['Failed rifts' if failed else 'Rifts'].append(_feature('ContinentalRift', geometry, '%s into plate %d' % (r['kind'].lower(), r['new']),
                                         r['plate'], begin, end,
                                         'split plate %d at %.0f Ma' % (r['new'], history.project.age_ma(r['time_myr']))))
        for belt in getattr(snapshot, 'orogens', []):
            for polygon in _polygons(belt['shape'], np.eye(3)):
                out['%s orogenies' % belt['kind'].capitalize()].append(_feature(
                    'OrogenicBelt', polygon, '%s orogeny' % belt['kind'], belt['plate'], begin, end,
                    'since %.0f Ma' % history.project.age_ma(belt['since_myr'])))
        for h, plate, points, first, last in plumes.trails(getattr(snapshot, 'hotspots', [])):
            geometry = pygplates.PolylineOnSphere(_points_on_sphere(points))
            out['Hotspot trails'].append(_feature('UnclassifiedFeature', geometry, '%s trail' % h['name'], plate, begin,
                                                  end, 'from %.0f to %.0f Ma' % (history.project.age_ma(first),
                                                                                  history.project.age_ma(last))))
        for f in getattr(snapshot, 'features', []):
            if f.kind == 'lip':
                young = snapshot.time_myr - f.born_myr <= plumes.ACTIVE_LIP_MYR + 1e-6
                for polygon in _polygons(f.shape, np.eye(3)):
                    out['Active LIPs' if young else 'Former LIPs'].append(_feature(
                        'LargeIgneousProvince', polygon, f.origin, f.plate, begin, end, _story(f, history.project)))
                continue
            layer, kind = FEATURE_LAYERS[f.kind]
            for polygon in _polygons(f.shape, np.eye(3)):
                out[layer].append(_feature(kind, polygon, '%s %d' % (f.kind.replace('_', ' '), f.ident), f.plate,
                                           begin, end, _story(f, history.project)))
        for b in world.boundaries:
            for kind, plate, line in _runs(b):
                feature_type, layer = BOUNDARY_TYPES.get(kind, ('UnclassifiedFeature', 'Transform'))
                back = world.frame(plate).T
                geometry = pygplates.PolylineOnSphere(_points_on_sphere(sphere.rotate(line, back)))
                feature = _feature(feature_type, geometry, '%s %d-%d' % (kind, b.left, b.right),
                                   plate, begin, end)
                if kind in ('ridge', 'transform'):
                    feature.set_left_plate(b.left)
                    feature.set_right_plate(b.right)
                    if kind == 'ridge':
                        feature.set_conjugate_plate_id(b.right)
                if kind == 'trench':
                    # The overriding plate carries the trench; its polarity says which side it is on.
                    feature.set_enumeration(pygplates.PropertyName.gpml_subduction_polarity,
                                            'Left' if plate == b.left else 'Right')
                out[layer].append(feature)
    _hotspots(history, out)
    return {name: _merged(features) for name, features in out.items()}


FEATURE_LAYERS = {
    'island_arc': ('Island arcs', 'IslandArc'),
    'continental_arc': ('Continental arcs', 'OrogenicBelt'),
    'terrane': ('Accreted terranes', 'ContinentalCrust'),
}


def _story(feature, project) -> str:
    """Where a feature came from and where it went, in ages, for its GPlates description."""
    text = '%s, from %.0f Ma' % (feature.origin, project.age_ma(feature.born_myr))
    for a in feature.accreted:
        text += '; accreted from plate %d to plate %d at %.0f Ma' % (a['from_plate'], a['to_plate'],
                                                                     project.age_ma(a['time_myr']))
    return text


def _runs(b):
    """(kind, carrying plate, points) for each run of one kind (and one overriding plate) along `b`."""
    runs = []
    n = len(b.kinds)
    start = 0
    for k in range(1, n + 1):
        if k == n or b.kinds[k] != b.kinds[start] or b.overriding[k] != b.overriding[start]:
            kind, over = b.kinds[start], b.overriding[start]
            plate = over if kind == 'trench' and over is not None else b.left
            runs.append((kind, plate, b.points[start:k + 1]))
            start = k
    return runs


def export(history, directory) -> list:
    """Write the rotation file and every layer into `directory`; returns the files written."""
    os.makedirs(directory, exist_ok=True)
    written = []
    rotation = os.path.join(directory, 'rotation.rot')
    write_rotations(history, rotation)
    written.append(rotation)
    for name, features in layers(history).items():
        path = os.path.join(directory, name + '.gpml')
        pygplates.FeatureCollection(features).write(path)
        written.append(path)
    return written
