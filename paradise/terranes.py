"""Island arcs that grow behind trenches, and terranes that accrete instead of being subducted.

Both are bookkeeping rules, not mechanics:

- **Arcs.** Behind every trench that has lasted a whole step and closes at least
  ARC_MIN_CLOSING_CM_YR, a band ARC_HALF_WIDTH_KM either side of a line ARC_OFFSET_KM inboard of the
  trench (on the overriding plate) is volcanic arc. Where that band lies on the overriding plate's
  ocean it becomes new, thin continental crust: an island arc. Where it lies on continent it is a
  continental (Andean) arc, recorded but adding no crust. An arc left behind when its trench moves
  away stays where it is, as a remnant arc; arcs do not ride their trenches across oceans.
- **Terranes.** Before each sub-step, any piece of continental crust no larger than
  TERRANE_MAX_KM2 (an island arc, a microcontinent) that will reach a trench its plate goes down
  is taken off the downgoing plate, with the ocean around it to the trench, and given to the
  overriding plate. The boundaries are derived again, keeping their state by position, so the
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
ARC_HALF_WIDTH_KM = 30.
ARC_THICKNESS_KM = 25.       # juvenile arc crust, thinner than mature continent
ARC_GROWTH_KM_PER_MYR = .2   # an arc still fed by its trench thickens towards mature crust
ARC_MATURE_KM = 35.
ARC_MIN_CLOSING_CM_YR = 1.
ARC_MIN_KM2 = 1000.          # an arc smaller than this is not worth keeping
TERRANE_MAX_KM2 = 1e6        # about Japan and Kamchatka together; anything larger collides
ACCRETE_MARGIN_KM = 50.      # how far ahead of a sub-step's motion a terrane is caught
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
    normal = sphere.unit(np.vstack((segment[:1], segment[:-1] + segment[1:], segment[-1:])))
    angle = (km if to_left else -km) / sphere.RADIUS_KM
    return sphere.unit(points * np.cos(angle) + normal * np.sin(angle))


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
def grow_arcs(world, trenches: list, before: list, features: list, time_myr: float, step_myr: float = 10.) -> list:
    """Add this step's arcs to `world` and `features`; returns events for arcs that are new.

    `trenches` are the step's trench records (paradise.events.trenches) and `before` the previous
    step's: only a trench of the same plates and polarity in both has been going down long enough to
    feed an arc.
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
        band = sphere.intersection(sphere.strip(front, ARC_HALF_WIDTH_KM), world.regions[over])
        if sphere.area_km2(band) < ARC_MIN_KM2:
            continue
        frame = world.frame(over)
        local = sphere.rotate_geography(band, frame.T, tolerance=1e-4)
        new_arc, on_continent = [], []
        kept = []
        for c in world.crust:
            if c.plate != over:
                kept.append(c)
                continue
            inside = sphere.intersection(c.shape, local)
            if sphere.area_km2(inside) < sphere.DETAIL_KM2:
                kept.append(c)
                continue
            if c.kind == 'ocean':
                c.shape = sphere.difference(c.shape, inside)
                new_arc.append(inside)
            else:
                on_continent.append(inside)
                if c.thickness_km < ARC_MATURE_KM and _is_arc(features, over, c.shape):
                    # Arc crust still fed by its trench matures; its thickness mesh is built again
                    # from the new mean.
                    c.thickness_km = min(ARC_MATURE_KM, c.thickness_km + ARC_GROWTH_KM_PER_MYR * step_myr)
                    c.mesh = None
            kept.append(c)
        world.crust = [c for c in kept if sphere.area_km2(c.shape) >= sphere.DETAIL_KM2]
        pair = '%d and %d' % tuple(t['plates'])
        if new_arc:
            shape = sphere.union(*new_arc)
            if sphere.area_km2(shape) >= ARC_MIN_KM2:
                world.crust.append(Crust(over, 'continent', shape, None, ARC_THICKNESS_KM))
                if _merge_into(features, 'island_arc', over, shape) is None:
                    ident = _next_ident(features)
                    features.append(Feature(ident, 'island_arc', over, shape, time_myr,
                                            'island arc of the trench between plates %s' % pair))
                    events.append(dict(kind='arc_born', plate=over, plates=list(t['plates']), feature=ident,
                                       area_km2=round(sphere.area_km2(shape))))
        if on_continent:
            # Arcs already there (island arcs, terranes) are not a new continental arc.
            shape = sphere.union(*on_continent)
            for f in features:
                if f.plate == over and f.kind in ('island_arc', 'terrane'):
                    shape = sphere.difference(shape, f.shape)
            if sphere.area_km2(shape) >= ARC_MIN_KM2 and _merge_into(features, 'continental_arc', over, shape) is None:
                features.append(Feature(_next_ident(features), 'continental_arc', over, shape, time_myr,
                                        'continental arc of the trench between plates %s' % pair))
    return events


def _is_arc(features, plate, shape) -> bool:
    """Whether most of a crust piece is an island arc of this plate."""
    area = sphere.area_km2(shape)
    arcs = [f.shape for f in features if f.plate == plate and f.kind == 'island_arc']
    return bool(arcs) and area > 0. and sphere.area_km2(sphere.intersection(shape, sphere.union(*arcs))) > .5 * area


def _merge_into(features, kind, plate, shape):
    """Join `shape` to a feature of this kind on this plate that it touches; the feature, or None."""
    for f in features:
        if f.kind == kind and f.plate == plate and _distance_km(f.shape, shape) < 1.:
            f.shape = sphere.union(f.shape, shape)
            return f
    return None


# ----------------------------------------------------------------------------------------- terranes
def accrete(world, features: list, rotations: dict, dt: float) -> list:
    """Hand terranes about to reach a trench to the overriding plate; returns 'accreted' events."""
    events = []
    for _ in range(8):   # each accretion changes the boundaries, so look again until none is left
        found = _next_accretion(world, rotations, dt)
        if found is None:
            break
        events.append(_transfer(world, features, rotations, *found))
    return events


def _next_accretion(world, rotations, dt):
    barred = getattr(world, 'no_accretion', ())
    for b, start, stop, over, down in _trench_runs(world):
        if down not in world.regions or over not in world.regions or frozenset((over, down)) in barred:
            continue
        opening, _ = segment_motion(b.points[start:stop + 1], rotations.get(b.left, np.zeros(3)),
                                    rotations.get(b.right, np.zeros(3)))
        closing = float(np.max(-opening, initial=0.))
        if closing <= 0.:
            continue
        reach = closing * dt + ACCRETE_MARGIN_KM
        line = _line(b.points[start:stop + 1])
        frame = world.frame(down)
        for c in world.crust:
            if c.plate != down or c.kind != 'continent':
                continue
            try:
                parts, _ = sphere.parts(sphere.rotate_geography(c.shape, frame, tolerance=1e-4))
            except (ValueError, RuntimeError):
                continue
            largest = getattr(world, 'terrane_max_km2', TERRANE_MAX_KM2)   # a project's own limit
            for part in parts:
                area = sphere.area_km2(part)
                if area <= largest and area >= sphere.DETAIL_KM2:
                    distance = _distance_km(part, line)
                    if distance <= reach:
                        others = [p for p in parts if p is not part]
                        return part, distance, over, down, others
    return None


def _transfer(world, features, rotations, part, distance, over, down, others) -> dict:
    """Give `part` (world frame), and the ocean round it out to the trench, from `down` to `over`."""
    rings = sphere.polygons(part)
    outline = np.vstack([np.vstack((r, r[:1])) for poly in rings for r in poly[:1]])
    zone = sphere.union(part, sphere.strip(outline, distance + ACCRETE_MARGIN_KM))
    zone = sphere.intersection(zone, world.regions[down])
    for other in others:   # never a bite out of the plate's other continents
        zone = sphere.difference(zone, other)
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
            moved.append(Crust(over, c.kind, sphere.rotate_geography(world_taken, to_over, tolerance=1e-4),
                               c.born_myr, c.thickness_km))
        kept.append(c)
    world.crust = [c for c in kept if sphere.area_km2(c.shape) >= sphere.DETAIL_KM2] + moved
    area = sphere.area_km2(part)
    # The features on the terrane go with it and become (or stay) a terrane.
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
            carried.append(f.ident)
    if not carried:
        ident = _next_ident(features)
        features.append(Feature(ident, 'terrane', over, sphere.rotate_geography(part, to_over, tolerance=1e-4),
                                world.time_myr, 'microcontinent from plate %d' % down, [dict(record)]))
        carried.append(ident)
    world.regions = world._from_crust()
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
        # Every feature lies on continental crust; whatever of it is no longer there has gone.
        continent = [c.shape for c in world.crust if c.plate == f.plate and c.kind == 'continent']
        f.shape = sphere.intersection(f.shape, sphere.union(*continent)) if continent else sphere.EMPTY
        if f.area_km2() >= ARC_MIN_KM2:
            result.append(f)
    return result
