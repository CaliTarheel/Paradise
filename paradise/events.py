"""Events a person declares at the start of a step, and the report on the trenches they make.

A declaration names two plates and a stretch of the boundary between them, as lon/lat points
picked on the map at the time of the step. It applies to every segment of that boundary within
REACH_KM of the stretch, so replaying the log finds the same segments: the world at the start of
the step is the same every time.

- `subduction`: the segments become a trench with the named plate overriding. The core keeps a
  trench a trench as the boundary moves (the overriding plate is carried with each segment), so
  the declaration lasts until it is ended or a continent arrives and the trench becomes a collision.
  Declaring it again with the other plate over it flips its polarity.
- `end_subduction`: the segments stop being a trench and are typed from the motion again, with
  their closing so far forgotten.

Two more declarations split a plate (paradise.rifting), and may give the new plate a pole:

- `rift`: a line drawn across a plate; the side named (`new_side`, left or right walking the line
  as drawn) becomes a new plate.
- `backarc`: behind a stretch of trench, the band of the overriding plate out to the arc and a
  little beyond becomes a new plate, which can roll back towards the trench.

And three decide whether two plates that collide become one (the core welds them once their
collision has closed SUTURE_KM):

- `weld`: weld them now, whatever their boundary: every segment between them becomes a collision
  closed SUTURE_KM, and the core welds them before anything moves.
- `keep_apart`: never weld them, from this step on; they go on colliding.
- `allow_weld`: undo `keep_apart`.
- `no_accretion`: continent reaching a trench between these plates collides, however small, from
  this step on (paradise.terranes reads the world's `no_accretion` set); `allow_accretion` undoes it.

Plumes (paradise.plumes): `hotspot` at a point, `end_hotspot` by number, and `lip` at a point
with a radius, leaving a hotspot behind unless `with_hotspot` is false.
"""
from __future__ import annotations

import numpy as np

from paradise.core import sphere
from paradise.core.build import MINIMUM_RUN, densify
from paradise.core.regions import SUTURE_KM
from paradise.core.topology import segment_motion

REACH_KM = 60.
KINDS = ('subduction', 'end_subduction', 'rift', 'backarc', 'weld', 'keep_apart', 'allow_weld',
         'no_accretion', 'allow_accretion',
         'hotspot', 'end_hotspot', 'lip')
PLUMES = ('hotspot', 'end_hotspot', 'lip')   # paradise.plumes, applied by the history
PAIRS = ('weld', 'keep_apart', 'allow_weld', 'no_accretion', 'allow_accretion')   # about two plates, not a stretch
SPLITS = ('rift', 'backarc')
# A plate carried by a slab moves fast (Scotese's first rule of thumb): below this at its trench,
# the downgoing plate's motion is worth a second look.
SLOW_SLAB_CM_YR = 2.
LONG_TRENCH_KM = 1000.


def validate(event: dict) -> dict:
    """The event in a canonical form, or ValueError saying what is wrong with it."""
    kind = event.get('kind')
    if kind not in KINDS:
        raise ValueError('unknown event %r' % kind)
    pole = event.get('pole')
    if pole is not None:
        pole = [float(v) for v in pole]
        if len(pole) != 3 or not -90. <= pole[0] <= 90.:
            raise ValueError('a pole is latitude, longitude and a rate in degrees per Myr')
    if kind in ('hotspot', 'lip'):
        lon, lat = (float(v) for v in event.get('point', ()))
        if not -90. <= lat <= 90.:
            raise ValueError('a %s needs a point on the map' % kind)
        result = dict(kind=kind, point=[lon, lat])
        if event.get('name'):
            result['name'] = str(event['name'])
        if kind == 'lip':
            radius = float(event.get('radius_km', 800.))
            if not 50. <= radius <= 3000.:
                raise ValueError('a LIP is between 50 and 3000 km in radius')
            result['radius_km'] = radius
            result['with_hotspot'] = bool(event.get('with_hotspot', True))
        return result
    if kind == 'end_hotspot':
        return dict(kind=kind, hotspot=int(event['hotspot']))
    if kind == 'rift':
        points = [[float(lon), float(lat)] for lon, lat in event.get('points', [])]
        if len(points) < 2:
            raise ValueError('a rift needs a line of two points or more')
        side = event.get('new_side', 'left')
        if side not in ('left', 'right'):
            raise ValueError('the new plate of a rift is on its left or its right')
        result = dict(kind=kind, plate=int(event['plate']), points=points, new_side=side)
        if pole is not None:
            result['pole'] = pole
        return result
    plates = [int(p) for p in event.get('plates', [])]
    if len(plates) != 2 or plates[0] == plates[1]:
        raise ValueError('a %s names two different plates' % kind.replace('_', ' '))
    if kind in PAIRS:
        return dict(kind=kind, plates=sorted(plates))
    points = [[float(lon), float(lat)] for lon, lat in event.get('points', [])]
    if len(points) < 2:
        raise ValueError('a %s needs a stretch of boundary (two points or more)' % kind.replace('_', ' '))
    result = dict(kind=kind, plates=sorted(plates), points=points)
    if pole is not None and kind == 'backarc':
        result['pole'] = pole
    if kind in ('subduction', 'backarc'):
        over = int(event.get('overriding', -1))
        if over not in plates:
            raise ValueError('the overriding plate must be one of plates %d and %d' % tuple(plates))
        result['overriding'] = over
    return result


def apply(world, events: list) -> list:
    """Apply declarations to `world`'s boundaries; returns notes for any that found nothing."""
    notes = []
    for event in events:
        event = validate(event)
        if event['kind'] in SPLITS or event['kind'] in PLUMES:
            continue   # paradise.rifting and paradise.plumes, applied by the history
        if event['kind'] in PAIRS:
            notes += _pair(world, event)
            continue
        left, right = event['plates']
        line = sphere.from_lonlat(*np.array(event['points']).T)
        dense = densify(line, spacing_km=REACH_KM / 3.)
        matched = 0
        for b in world.boundaries:
            if (b.left, b.right) != (left, right) or len(b.points) < 2:
                continue
            mids = b.mids()
            near = np.min(sphere.arc_km(mids[:, None, :], dense[None, :, :]), axis=1) <= REACH_KM
            for k in np.flatnonzero(near):
                if event['kind'] == 'subduction':
                    b.kinds[k], b.overriding[k] = 'trench', event['overriding']
                else:
                    b.kinds[k], b.overriding[k] = 'converging', None
                b.closed_km[k] = 0.
            matched += int(near.sum())
        if 0 < matched < MINIMUM_RUN:
            notes.append(dict(kind='declaration_short', plates=[left, right], text=(
                'The declared %s between plates %d and %d covers only %d boundary segment%s; runs shorter '
                'than %d take the type of their neighbours. Pick a longer stretch.'
                % (event['kind'].replace('_', ' '), left, right, matched, '' if matched == 1 else 's', MINIMUM_RUN))))
        if not matched:
            notes.append(dict(kind='declaration_missed', plates=[left, right],
                              text='The declared %s between plates %d and %d found no boundary there any more.'
                                   % (event['kind'].replace('_', ' '), left, right)))
    return notes


def _pair(world, event) -> list:
    pair = frozenset(event['plates'])
    if event['kind'] in ('no_accretion', 'allow_accretion'):
        barred = set(getattr(world, 'no_accretion', ()))
        world.no_accretion = barred | {pair} if event['kind'] == 'no_accretion' else barred - {pair}
        return []
    kept = set(getattr(world, 'kept_apart', ()))
    if event['kind'] == 'keep_apart':
        world.kept_apart = kept | {pair}
        return []
    world.kept_apart = kept - {pair}
    if event['kind'] == 'allow_weld':
        return []
    left, right = event['plates']
    found = False
    for b in world.boundaries:
        if (b.left, b.right) == (left, right):
            b.kinds = ['collision'] * len(b.kinds)
            b.overriding = [None] * len(b.kinds)
            b.closed_km = np.full(len(b.kinds), SUTURE_KM)
            found = found or len(b.kinds) > 0
    if not found:
        return [dict(kind='declaration_missed', plates=[left, right],
                     text='Plates %d and %d do not touch, so they cannot be welded.' % (left, right))]
    return []


def collisions(world) -> list:
    """Pairs of plates colliding now: [{plates, length_km, kept_apart}]."""
    kept = getattr(world, 'kept_apart', set())
    found = {}
    for b in world.boundaries:
        if len(b.points) < 2:
            continue
        lengths = sphere.arc_km(b.points[:-1], b.points[1:])
        km = float(sum(l for k, l in zip(b.kinds, lengths) if k == 'collision'))
        if km > 0.:
            found[(b.left, b.right)] = found.get((b.left, b.right), 0.) + km
    return [dict(plates=list(pair), length_km=round(km), kept_apart=frozenset(pair) in kept)
            for pair, km in sorted(found.items())]


def trenches(world, vectors: dict, step_myr: float) -> list:
    """Every trench as it stands, measured with the poles of the step that led here.

    One record per run of trench segments with one overriding plate: its plates, length, mean
    closing rate (negative when the plates pull apart across it), the age of the crust going down
    and the area consumed over a step at that rate.
    """
    ages = world._slab_ages()
    records = []
    for b in world.boundaries:
        if len(b.points) < 2:
            continue
        omega_left = np.asarray(vectors.get(b.left, np.zeros(3)), float)
        omega_right = np.asarray(vectors.get(b.right, np.zeros(3)), float)
        opening, _ = segment_motion(b.points, omega_left, omega_right)
        lengths = sphere.arc_km(b.points[:-1], b.points[1:])
        start = 0
        n = len(b.kinds)
        for k in range(1, n + 1):
            if k < n and (b.kinds[k], b.overriding[k]) == (b.kinds[start], b.overriding[start]):
                continue
            if b.kinds[start] == 'trench' and b.overriding[start] in (b.left, b.right):
                run = slice(start, k)
                length = float(lengths[run].sum())
                if length > 0.:
                    closing = float(-(opening[run] * lengths[run]).sum() / length)   # km/Myr
                    slab = [ages[(id(b), j)] for j in range(start, k) if (id(b), j) in ages]
                    over = b.overriding[start]
                    down = b.right if over == b.left else b.left
                    middle = (start + k) // 2
                    sides = world._sides(b, middle)
                    down_crust = sides[0] if down == b.left else sides[1]
                    records.append(dict(
                        plates=[b.left, b.right], overriding=over, downgoing=down,
                        length_km=round(length), closing_cm_yr=round(closing / 10., 2),
                        slab_age_myr=round(float(np.mean(slab)), 1) if slab else None,
                        consumed_km2_per_step=round(max(closing, 0.) * length * step_myr),
                        downgoing_crust=down_crust,
                        points=[[round(float(lon), 3), round(float(lat), 3)]
                                for lon, lat in zip(*sphere.to_lonlat(b.points[start:k + 1]))]))
            start = k
    return records


def trench_warnings(records: list, vectors: dict) -> list:
    """What a person might want to look at again in the trenches of a step."""
    warnings = []
    for t in records:
        if t['closing_cm_yr'] < 0.:
            warnings.append(dict(kind='trench_opening', plates=t['plates'], text=(
                'Plate %d pulls away from its trench with plate %d at %.1f cm/yr. The gap fills with new '
                'ocean split between the two plates; there is no back-arc spreading yet.'
                % (t['downgoing'], t['overriding'], -t['closing_cm_yr']))))
        if t['downgoing_crust'] == 'continent':
            warnings.append(dict(kind='continent_at_trench', plates=t['plates'], text=(
                'Continent on plate %d is at its trench under plate %d; once enough arrives the trench '
                'becomes a collision.' % (t['downgoing'], t['overriding']))))
        if t['length_km'] >= LONG_TRENCH_KM:
            middle = t['points'][len(t['points']) // 2]
            point = sphere.from_lonlat(*middle)
            omega = np.asarray(vectors.get(t['downgoing'], np.zeros(3)), float)
            speed = float(np.linalg.norm(np.cross(omega, point))) * sphere.RADIUS_KM / 10.
            if speed < SLOW_SLAB_CM_YR:
                warnings.append(dict(kind='slow_slab', plates=t['plates'], text=(
                    'Plate %d has %d km of trench pulling on it but moves at only %.1f cm/yr there; plates '
                    'carried by a slab usually move fast.' % (t['downgoing'], t['length_km'], speed))))
    return warnings
