"""Island arcs that grow behind trenches, and terranes that accrete instead of being subducted.

Both are bookkeeping rules, not mechanics:

- **Arcs.** Behind every trench that has lasted a whole step and closes at least
  ARC_MIN_CLOSING_CM_YR, the volcanic front ARC_OFFSET_KM inboard of the trench (on the overriding
  plate) is cut into stretches, and each stretch with no arc yet gets one: a band ARC_HALF_WIDTH_KM
  either side of it, an island arc over ocean or a continental (Andean) arc over continent. Arcs
  are features, as in the hand-built Artifexia project, not crust: making them crust cut the
  ocean into hundreds of slivers. An arc left behind when its trench moves away stays where it is,
  as a remnant arc; arcs do not ride their trenches across oceans.
- **Terranes.** Before each sub-step, any piece of continental crust no larger than
  TERRANE_MAX_KM2 (an island arc, a microcontinent) that will reach a trench its plate goes down
  is taken off the downgoing plate, with the ocean around it to the trench, and given to the
  overriding plate. That ocean is welded on as continental crust, and the terrane is all of it. The boundaries are derived again, keeping their state by position, so the
  trench steps outboard of the new terrane and carries on. Larger continents still collide.

The arcs and terranes are kept as features beside the world, each in its plate's frame like the
cratons, so they keep their identity and history (where they were born, when and where they
accreted) while the core's crust pieces merge and split.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import spherely

from paradise.core import sphere
from paradise.core.topology import segment_motion
from paradise.core.world import Crust

ARC_OFFSET_KM = 150.         # trench to volcanic front, as in most arcs (about 100 to 200 km)
ARC_HALF_WIDTH_KM = 15.      # 30 km wide, as the hand-built Artifexia arcs are (median 28 km)
ARC_SEGMENT_KM, ARC_GAP_KM = 800., 150.   # Artifexia's arcs are separate, typically 800 km long
ARC_NEAR_KM = 100.           # a stretch of front with an arc this close has one already
ARC_MIN_CLOSING_CM_YR = 1.
ARC_MIN_KM2 = 1000.          # an arc smaller than this is not worth keeping
TERRANE_MAX_KM2 = 1e6        # about Japan and Kamchatka together; anything larger collides
# Crumbs smaller than this are not terranes: the core takes continent arriving at a trench in pieces
# this small as shortening, without collision.
TERRANE_MIN_KM2 = 15000.     # Artifexia's arcs are about 24,000 km2; its smallest terranes 60,000
ACCRETIONS_PER_SUBSTEP = 3   # more than this a sub-step were crumbs passed back and forth
# A terrane that has just accreted does not accrete again this soon. Where two plates meet along
# trenches of both polarities, it would otherwise pass from one plate to the other and straight back.
SETTLE_MYR = 10.
ACCRETE_MARGIN_KM = 50.      # how far ahead of a sub-step's motion a terrane is caught
# The ocean that accretes with a terrane (forearc, accretionary wedge, slivers of ocean floor) is
# welded on as continental crust this thick, and is part of the terrane: in the hand-built
# Artifexia project a terrane is a median eight times the area of the arcs it grew from.
ACCRETED_KM = 30.
STEPPED_KM = 5.              # boundary segments this close to an accreted zone are its new trench
_KM_PER_S2_METRE = sphere.RADIUS_KM / spherely.EARTH_RADIUS_METERS


@dataclass
class Feature:
    """An island arc, a continental arc or a terrane, in its plate's frame."""
    ident: int
    kind: str                # 'island_arc', 'continental_arc' or 'terrane'
    plate: int
    shape: object
    born_myr: float
    origin: str
    accreted: list = field(default_factory=list)   # [dict(time_myr, from_plate, to_plate, area_km2)]

    def area_km2(self) -> float:
        return sphere.area_km2(self.shape)


def _next_ident(features) -> int:
    return max((f.ident for f in features), default=0) + 1


def _distance_km(a, b) -> float:
    return float(spherely.distance(a, b)) * _KM_PER_S2_METRE


def _line(points):
    lon, lat = sphere.to_lonlat(np.asarray(points, float))
    return spherely.create_linestring(list(zip(lon.tolist(), lat.tolist())))


def offset_line(points, km, to_left: bool) -> np.ndarray:
    """`points` moved `km` to one side, square to the line at each vertex."""
    points = np.asarray(points, float)
    # Repeated points have no direction between them.
    points = points[np.r_[True, sphere.arc_km(points[:-1], points[1:]) > 1e-3]]
    if len(points) < 2:
        return points
    segment = sphere.unit(np.cross(points[:-1], points[1:]))   # left of each segment, seen from outside
    summed = np.vstack((segment[:1], segment[:-1] + segment[1:], segment[-1:]))
    # Where the line turns straight back, its two segments' normals cancel: take the first one's.
    back = np.linalg.norm(summed, axis=1) < 1e-6
    summed[back] = np.vstack((segment[:1], segment))[back]
    normal = sphere.unit(summed)
    angle = (km if to_left else -km) / sphere.RADIUS_KM
    moved = sphere.unit(points * np.cos(angle) + normal * np.sin(angle))
    return moved[np.all(np.isfinite(moved), axis=1)]


def _trench_runs(world):
    """(boundary, start, stop, overriding, downgoing) for every run of trench with one overriding plate."""
    for b in world.boundaries:
        n = len(b.kinds)
        start = 0
        for k in range(1, n + 1):
            if k < n and (b.kinds[k], b.overriding[k]) == (b.kinds[start], b.overriding[start]):
                continue
            over = b.overriding[start]
            if b.kinds[start] == 'trench' and over in (b.left, b.right):
                yield b, start, k, over, (b.right if over == b.left else b.left)
            start = k


# --------------------------------------------------------------------------------------------- arcs
def _chunks(points, length_km, gap_km):
    """`points` cut into stretches about `length_km` long with `gap_km` between them."""
    lengths = np.r_[0., np.cumsum(sphere.arc_km(points[:-1], points[1:]))]
    total = lengths[-1]
    if total <= length_km:
        return [points]
    count = max(1, int(round((total + gap_km) / (length_km + gap_km))))
    stride = (total + gap_km) / count
    chunks = []
    for k in range(count):
        lo, hi = k * stride, k * stride + stride - gap_km
        inside = (lengths >= lo) & (lengths <= hi)
        if inside.sum() >= 2:
            chunks.append(points[inside])
    return chunks


def grow_arcs(world, trenches: list, before: list, features: list, time_myr: float, step_myr: float = 10.) -> list:
    """Add this step's new arcs to `features`; returns events for them.

    `trenches` are the step's trench records (paradise.events.trenches) and `before` the previous
    step's: only a trench of the same plates and polarity in both has been going down long enough to
    feed an arc. Its volcanic front is cut into stretches ARC_SEGMENT_KM long, as the hand-built
    Artifexia arcs are, and each stretch with no arc of its plate within ARC_NEAR_KM gets one. Arcs
    are features on the plate, not crust: they leave the core's geometry alone.
    """
    lasting = {(tuple(t['plates']), t['overriding']) for t in before}
    events = []
    for t in trenches:
        over = t['overriding']
        if (tuple(t['plates']), over) not in lasting or t['closing_cm_yr'] < ARC_MIN_CLOSING_CM_YR:
            continue
        if over not in world.regions or len(t['points']) < 2:
            continue
        points = sphere.from_lonlat(*np.array(t['points']).T)
        front = offset_line(points, ARC_OFFSET_KM, to_left=(over == t['plates'][0]))
        if len(front) < 2:
            continue
        frame = world.frame(over)
        continent = [c.shape for c in world.crust if c.plate == over and c.kind == 'continent']
        continent = sphere.union(*continent) if continent else None
        pair = '%d and %d' % tuple(t['plates'])
        for chunk in _chunks(front, ARC_SEGMENT_KM, ARC_GAP_KM):
            band = sphere.intersection(sphere.strip(chunk, ARC_HALF_WIDTH_KM), world.regions[over])
            area = sphere.area_km2(band)
            if area < ARC_MIN_KM2:
                continue
            local = sphere.rotate_geography(band, frame.T, tolerance=1e-4)
            on_land = sphere.area_km2(sphere.intersection(local, continent)) if continent is not None else 0.
            kind = 'continental_arc' if on_land > .5 * area else 'island_arc'
            if any(f.plate == over and f.kind == kind and _distance_km(f.shape, local) < ARC_NEAR_KM for f in features):
                continue   # an arc is there already
            ident = _next_ident(features)
            features.append(Feature(ident, kind, over, local, time_myr,
                                    '%s of the trench between plates %s' % (kind.replace('_', ' '), pair)))
            if kind == 'island_arc':
                events.append(dict(kind='arc_born', plate=over, plates=list(t['plates']), feature=ident,
                                   area_km2=round(area)))
    return events


# ----------------------------------------------------------------------------------------- terranes
def accrete(world, features: list, rotations: dict, dt: float) -> list:
    """Hand terranes about to reach a trench to the overriding plate; returns 'accreted' events."""
    import copy
    events = []
    for _ in range(ACCRETIONS_PER_SUBSTEP):   # each changes the boundaries, so look again each time
        found = _next_accretion(world, features, rotations, dt)
        if found is None:
            break
        # All or nothing: rebuilding the plates after a transfer can trip the core's checks (a
        # sliver where the accreted crust meets its new plate), and a half-made transfer would
        # spoil the whole step. Undone, the terrane simply waits, or goes down with its plate.
        saved = (copy.deepcopy(world.crust), dict(world.regions), copy.deepcopy(world.boundaries),
                 copy.deepcopy(features))
        try:
            events.append(_transfer(world, features, rotations, *found))
        except (ValueError, RuntimeError) as e:
            world.crust, world.regions, world.boundaries = saved[0], saved[1], saved[2]
            features[:] = saved[3]
            events.append(dict(kind='accretion_failed', plate=found[3], into=found[2], text=str(e)))
            break
    return events


def _settling_shapes(world, features) -> list:
    """Terranes that accreted less than SETTLE_MYR ago, in the world frame."""
    return [sphere.rotate_geography(f.shape, world.frame(f.plate), tolerance=1e-4) for f in features
            if f.kind == 'terrane' and f.accreted and world.time_myr - f.accreted[-1]['time_myr'] < SETTLE_MYR
            and f.plate in world.regions]


def _settling(settling, part) -> bool:
    """Whether `part` (world frame) is mostly one of the `settling` terranes."""
    area = sphere.area_km2(part)
    return any(sphere.area_km2(sphere.intersection(shape, part)) > .5 * area for shape in settling)


def _next_accretion(world, features, rotations, dt):
    barred = getattr(world, 'no_accretion', ())
    largest = getattr(world, 'terrane_max_km2', TERRANE_MAX_KM2)   # a project's own limit
    pieces = {}   # downgoing plate -> its continental pieces in the world frame, worked out once
    settling = _settling_shapes(world, features)

    def continental_parts(plate):
        if plate not in pieces:
            frame = world.frame(plate)
            found = []
            for c in world.crust:
                if c.plate == plate and c.kind == 'continent':
                    try:
                        parts, _ = sphere.parts(sphere.rotate_geography(c.shape, frame, tolerance=1e-4))
                    except (ValueError, RuntimeError):
                        continue
                    found += parts
            pieces[plate] = found
        return pieces[plate]

    def arc_shapes(plate):
        # Island arcs are features on the ocean, not crust, and accrete like continental pieces.
        frame = world.frame(plate)
        return [sphere.rotate_geography(f.shape, frame, tolerance=1e-4) for f in features
                if f.plate == plate and f.kind == 'island_arc' and not sphere.is_empty(f.shape)]

    for b, start, stop, over, down in _trench_runs(world):
        if down not in world.regions or over not in world.regions or frozenset((over, down)) in barred:
            continue
        opening, _ = segment_motion(b.points[start:stop + 1], rotations.get(b.left, np.zeros(3)),
                                    rotations.get(b.right, np.zeros(3)))
        closing = float(np.max(-opening, initial=0.))
        if closing <= 0.:
            continue
        parts = continental_parts(down) + arc_shapes(down)
        if not parts:
            continue
        reach = closing * dt + ACCRETE_MARGIN_KM
        line = _line(b.points[start:stop + 1])
        for part in parts:
            area = sphere.area_km2(part)
            if TERRANE_MIN_KM2 <= area <= largest:
                distance = _distance_km(part, line)
                if distance <= reach and not _settling(settling, part):
                    return part, distance, over, down, [p for p in parts if p is not part]
    return None


def _transfer(world, features, rotations, part, distance, over, down, others) -> dict:
    """Give `part` (world frame), and the ocean round it out to the trench, from `down` to `over`."""
    rings = sphere.polygons(part)
    outline = np.vstack([np.vstack((r, r[:1])) for poly in rings for r in poly[:1]])
    zone = sphere.union(part, sphere.strip(outline, distance + ACCRETE_MARGIN_KM))
    zone = sphere.intersection(zone, world.regions[down])
    for other in others:   # never a bite out of the plate's other continents
        zone = sphere.difference(zone, other)
    if sphere.area_km2(zone) < TERRANE_MIN_KM2:
        raise ValueError('nothing to accrete there')   # undone by accrete(), like any failed transfer
    to_down, to_over = world.frame(down).T, world.frame(over).T
    local = sphere.rotate_geography(zone, to_down, tolerance=1e-4)
    moved, kept = [], []
    for c in world.crust:
        if c.plate != down:
            kept.append(c)
            continue
        taken = sphere.intersection(c.shape, local)
        if sphere.area_km2(taken) >= sphere.DETAIL_KM2:
            c.shape = sphere.difference(c.shape, taken)
            # The piece's frame changes, so a continent's thickness mesh is built again from its mean.
            world_taken = sphere.rotate_geography(taken, world.frame(down), tolerance=1e-4)
            # Ocean taken with the terrane is welded on with it, as accreted crust.
            kind, born, thick = (c.kind, c.born_myr, c.thickness_km) if c.kind == 'continent' else \
                ('continent', None, ACCRETED_KM)
            moved.append(Crust(over, kind, sphere.rotate_geography(world_taken, to_over, tolerance=1e-4), born, thick))
        kept.append(c)
    world.crust = [c for c in kept if sphere.area_km2(c.shape) >= sphere.DETAIL_KM2] + moved
    area = sphere.area_km2(zone)
    zone_over = sphere.rotate_geography(zone, to_over, tolerance=1e-4)
    # The features on the terrane go with it and become (or stay) a terrane: all of the zone.
    record = dict(time_myr=world.time_myr, from_plate=down, to_plate=over, area_km2=round(area))
    carried = []
    for f in features:
        if f.plate != down:
            continue
        world_shape = sphere.rotate_geography(f.shape, world.frame(down), tolerance=1e-4)
        if sphere.area_km2(sphere.intersection(world_shape, zone)) >= .5 * f.area_km2():
            f.plate, f.shape = over, sphere.rotate_geography(world_shape, to_over, tolerance=1e-4)
            f.kind = 'terrane'
            f.accreted.append(dict(record))
            carried.append(f)
    if carried:
        # Arcs accreting together are one terrane (in Artifexia, a third of terranes grew from
        # two arcs or more), and it is the whole accreted zone.
        first = carried[0]
        first.shape = sphere.union(zone_over, *[f.shape for f in carried])
        if len(carried) > 1:
            first.origin = '%s, with %s' % (first.origin, ', '.join('%s %d' % (f.origin, f.ident) for f in carried[1:]))
        features[:] = [f for f in features if f not in carried[1:]]
        carried = [first.ident]
    else:
        ident = _next_ident(features)
        features.append(Feature(ident, 'terrane', over, zone_over, world.time_myr,
                                'microcontinent from plate %d' % down, [dict(record)]))
        carried = [ident]
    # Only the two plates changed: rebuilding every plate's region took half an accretion's time.
    for p in (down, over):
        shapes = [c.shape for c in world.crust if c.plate == p]
        if shapes:
            world.regions[p] = sphere.rotate_geography(sphere.union(*shapes), world.frame(p))
        else:
            world.regions.pop(p, None)
    world.boundaries = world._derive(world._index(), rotations)
    # The boundary round the new terrane can bulge further outboard than the core carries a
    # segment's state (INHERIT_KM), and would come back as an ordinary converging boundary. The
    # trench steps outboard with the terrane: every segment between the two plates along it is trench.
    pair = (min(down, over), max(down, over))
    for b in world.boundaries:
        if (b.left, b.right) != pair or len(b.points) < 2:
            continue
        lon, lat = sphere.to_lonlat(b.mids())
        mids = spherely.points(lon, lat)
        near = np.asarray(spherely.distance(mids, zone)) * _KM_PER_S2_METRE <= STEPPED_KM
        for k in np.flatnonzero(near):
            if b.kinds[k] != 'trench':
                b.kinds[k], b.overriding[k], b.closed_km[k] = 'trench', over, 0.
    return dict(kind='accreted', plates=[down, over], plate=down, into=over, features=carried, area_km2=round(area))


# ------------------------------------------------------------------------------------------ keeping
def tidy(world, features: list) -> list:
    """Features kept to their plate's continental crust, and moved with crust whose plate has gone."""
    plates = set(world.plates())
    result = []
    continents = {}   # plate -> union of its continental crust, worked out once (it was once a feature)

    def continent_of(plate):
        if plate not in continents:
            shapes = [c.shape for c in world.crust if c.plate == plate and c.kind == 'continent']
            continents[plate] = sphere.union(*shapes) if shapes else None
        return continents[plate]

    for f in features:
        if f.plate not in plates:
            # The plate was absorbed: the feature goes to the plate now under its middle.
            world_shape = sphere.rotate_geography(f.shape, world.frames.get(f.plate, np.eye(3)), tolerance=1e-4)
            c = spherely.centroid(world_shape)
            point = spherely.create_point(spherely.get_x(c), spherely.get_y(c))
            new = next((p for p in plates if spherely.contains(world.regions[p], point)), None)
            if new is None:
                continue
            f.plate, f.shape = new, sphere.rotate_geography(world_shape, world.frame(new).T, tolerance=1e-4)
        if f.kind == 'island_arc':
            # An island arc is a feature on its plate's ocean: whatever of it has left the plate
            # (gone down a trench, say) has gone.
            # Clipped in the world frame: the arc is small, its plate may be most of the sphere.
            frame = world.frame(f.plate)
            here = sphere.intersection(sphere.rotate_geography(f.shape, frame, tolerance=1e-4), world.regions[f.plate])
            f.shape = sphere.rotate_geography(here, frame.T, tolerance=1e-4) if not sphere.is_empty(here) else sphere.EMPTY
        else:
            # Every other feature lies on continental crust; whatever of it is no longer there has gone,
            # unless the crust went to another plate (the core passes pieces between plates), when the
            # feature goes with it.
            before = f.area_km2()
            continent = continent_of(f.plate)
            kept = sphere.intersection(f.shape, continent) if continent is not None else sphere.EMPTY
            if sphere.area_km2(kept) < .5 * before:
                world_shape = sphere.rotate_geography(f.shape, world.frame(f.plate), tolerance=1e-4)
                best, best_km2 = None, sphere.area_km2(kept)
                for p in plates:
                    if p == f.plate or continent_of(p) is None:
                        continue
                    here = sphere.intersection(sphere.rotate_geography(world_shape, world.frame(p).T, tolerance=1e-4),
                                               continent_of(p))
                    if sphere.area_km2(here) > best_km2:
                        best, best_km2 = (p, here), sphere.area_km2(here)
                if best is not None:
                    f.plate, kept = best
            f.shape = kept
        if f.area_km2() >= ARC_MIN_KM2:
            result.append(f)
    return result
