"""An advisor that proposes each step's poles and events from Scotese's rules of thumb.

It proposes; the history only records what is committed, as poles and declarations, exactly as if
a person had typed them. An autopilot (examples/thousand.py) commits its proposals step after step.

The rules:

- **Slab pull** (rules I and II). A plate with trenches where it goes down is pulled towards them:
  its rotation is the sum, over its downgoing trench segments, of the rotation that would carry the
  segment towards the overriding plate, weighted by length. Its speed grows with the share of its
  outline that is downgoing trench, from SLOW_CM_YR to SLAB_CM_YR, and continent on it slows it.
- **Ridge push** adds a weaker push away from the plate's ridges.
- **Undriven plates** keep one stage pole for STAGE_STEPS steps (Artifexia's stages are 50 Myr),
  more or less in the direction they were going. Oceanic ones barely move; continental ones drift at
  CONTINENTAL_CM_YR, which is what starts subduction in a Worldbuilding Pasta start.
- **Supercontinents break up** (rules III and XI). A plate carrying more than RIFT_CONTINENT_KM2 of
  continent that has not rifted for RIFT_TENURE_MYR may rift, along a great circle through the middle
  of its continent, preferring a suture if it has one, and the new plate moves away.
- **Continents override** (Worldbuilding Pasta's rule for its starting plates). Where a continental
  plate meets an oceanic one along a converging boundary or a trench, the advisor declares the
  trench with the continental plate overriding. The core's own polarity rule looks only 150 km
  either side, which on a continent's passive margin finds ocean on both, and so often sent the
  continental plate down, dragging its continent into the trench.
- **Large ocean plates split** along a new mid-ocean ridge (rule III: one ridge is all you need).
  An oceanic plate covering RIDGE_PLATE_SHARE of the sphere whose ocean is older than
  RIDGE_OCEAN_AGE_MYR on average rifts along a great circle through its point farthest from any
  continent, and the two parts move apart; slab pull then drives them into the trenches.
- **Old passive margins collapse.** A continental plate's ocean older than OLD_OCEAN_MYR, if there
  is OLD_OCEAN_KM2 of it, is detached as a new plate moving towards the continent, so that it goes
  down beneath it.
- Otherwise trenches form where plates converge (the project's automatic initiation), and a new
  pole for a plate is blended with its old one, so that motions change over a few steps.
"""
from __future__ import annotations

import math

import numpy as np
import spherely

from paradise.core import sphere
from paradise.history import pole_vector, vector_pole

SLOW_CM_YR, SLAB_CM_YR, MAX_CM_YR = 1., 9., 12.
# Continental plates are what moves first: in the Worldbuilding Pasta way of starting a world, one
# ocean plate and a plate per continent, nothing subducts until the continents push into the ocean.
CONTINENTAL_CM_YR = (3., 6.)
CONTINENTAL_SHARE = .3   # a plate this much continent or more is a continental plate
# A plate covering this much of the sphere or more (the one ocean plate of a Pasta start, or a
# supercontinent's) barely moves: every trench round it pulls a different way, and turning most of
# the sphere at slab speed would carry everything on it.
LARGE_SHARE, LARGE_CM_YR = .35, 1.
# A plate carrying this much continent drifts as a continental plate, however much ocean it has
# gathered at its ridges (the pieces of a rifted Pangaea grow large and mostly oceanic).
CONTINENTAL_KM2 = 10e6
POLE_STEER_DEG = 50.     # a drifting plate centred beyond this latitude is steered back
# A rifted plate goes on moving away from its parent for this long, at SEPARATION_CM_YR. With
# 60 Myr at 3 cm/yr, a Pangaea's pieces had parted by narrow seaways after 300 Myr.
SEPARATION_MYR, SEPARATION_CM_YR = 100., 4.
RIDGE_WEIGHT = .25
INERTIA = .5             # how much of a plate's old rotation each new one keeps
STAGE_STEPS = 5
RIFT_CONTINENT_KM2 = 10e6   # Artifexia's continents are mostly 10 to 50 million km2
RIFT_TENURE_MYR = 60.      # Artifexia's continents split every 50 to 100 Myr
RIFT_CHANCE = .35        # per step, once a plate may rift
RIFT_TRIES = 12          # directions tried for a rift that opens wherever it crosses the continent
RIDGE_CLEAR_KM = 300.     # a new ridge on a plate with a continent keeps this far from it
RIDGE_CAP_MIN_KM, RIDGE_CAP_MAX_KM = 800., 2500.   # the radius of ocean that comes away instead
RIFT_LAND_KM2 = 2e6       # land a rift is opened across (smaller pieces may close)
CIRCLE_STEP_DEG = .5     # the rift's great circle is looked at every 55 km
CROSSING_KM = 100.       # and crosses land wherever it comes this near it
# A large ocean plate splits along a new mid-ocean ridge, as a Pasta world adds ridges to its one
# ocean plate: without it, ocean leaves only where continents run over it, and grows ever older.
# Artifexia has 30 to 46 plates a stage: ocean plates are split far more often than one per globe.
RIDGE_PLATE_SHARE, RIDGE_OCEAN_AGE_MYR, RIDGE_TENURE_MYR = .08, 80., 40.
# Old passive margins collapse (Scotese IV: subduction must start somewhere). A continental plate's
# ocean is never subducted while it is the continent's own; a run at 600 Ma had a continental plate
# with half the sphere and ocean 300 Myr old on average. Ocean this old, this much of it, becomes a
# plate of its own, moving towards the continent, which then overrides it.
OLD_OCEAN_MYR, OLD_OCEAN_KM2, COLLAPSE_CM_YR = 160., 8e6, 4.


def _detached(event) -> bool:
    """Whether a rift event is a plate's old ocean detached (rifting.split_off names the new plate)."""
    return event.get('name', '').endswith('old ocean')


def _main_land(land):
    """The bodies of a continent a rift must open across: those of RIFT_LAND_KM2 or more.

    A supercontinent's plate gathers terranes and islands all round it; every great circle met one
    of them on its closing half, and the 127 Mkm2 supercontinent of the Pangaea run never rifted.
    A small piece on that half goes down the trench it meets, as such pieces do.
    """
    try:
        pieces, _ = sphere.parts(land)
    except ValueError:
        return land
    main = [p for p in pieces if sphere.area_km2(p) >= RIFT_LAND_KM2]
    return sphere.union(*main) if main else land


def _clear_of(land, centre, heading) -> float:
    """How near (km) the great circle through `centre` along `heading` comes to `land`."""
    around = np.radians(np.arange(0., 360., CIRCLE_STEP_DEG))
    points = np.outer(np.cos(around), centre) + np.outer(np.sin(around), heading)
    lon, lat = sphere.to_lonlat(points)
    near = np.asarray(spherely.distance(spherely.points(lon, lat), land))
    return float(near.min()) * sphere.RADIUS_KM / spherely.EARTH_RADIUS_METERS


def _opening_axis(land, centre, heading):
    """The axis that opens the great circle through `centre` along `heading` wherever it crosses land.

    A rotation opens a great circle on the half of it between its pole and the antipode, and closes
    it on the other half. The pole goes on the circle in the ocean, so that every crossing of the
    continent (`land`, a geography in the world frame) is on the opening half, and it turns the left
    side (the new plate) away. None when the crossings span more than half the circle.
    """
    side = np.cross(centre, heading)
    around = np.radians(np.arange(0., 360., CIRCLE_STEP_DEG))
    points = np.outer(np.cos(around), centre) + np.outer(np.sin(around), heading)
    lon, lat = sphere.to_lonlat(points)
    # Near land counts: a coast running along the circle, or a neck of land between two samples,
    # still meets it, and closing there welded the pieces again.
    near = np.asarray(spherely.distance(spherely.points(lon, lat), land)) * sphere.RADIUS_KM / spherely.EARTH_RADIUS_METERS
    on_land = near <= CROSSING_KM
    if not on_land.any():
        return None
    angles = np.sort(np.arctan2(points[on_land] @ heading, points[on_land] @ centre))
    gaps = np.diff(np.concatenate((angles, angles[:1] + 2. * math.pi)))
    k = int(np.argmax(gaps))
    if gaps[k] <= math.pi:
        return None
    # The crossings run from angles[k + 1] round to angles[k]: the pole goes before them, by half
    # what is left of the half circle, and its antipode as far after them.
    first = angles[(k + 1) % len(angles)]
    span = (angles[k] - first) % (2. * math.pi)
    t = first - .5 * (math.pi - span)
    pole = centre * math.cos(t) + heading * math.sin(t)
    middle = first + .5 * span
    probe = centre * math.cos(middle) + heading * math.sin(middle)
    return pole if np.dot(np.cross(pole, probe), side) > 0. else -pole


def _cm_yr_to_rad_myr(cm_yr, axis, centre):
    """The rotation rate (radians per Myr) about `axis` that moves `centre` at `cm_yr`.

    An axis near the centre would need a fast spin to move the centre at all; the rate is held to
    twice the plain one, so the plate's far edge never runs much faster than `cm_yr` either.
    """
    s = float(np.linalg.norm(np.cross(axis, centre)))
    return (cm_yr * 10. / sphere.RADIUS_KM) / max(s, .5)


def _centre(geography):
    c = spherely.centroid(geography)
    return sphere.from_lonlat(spherely.get_x(c), spherely.get_y(c))


class Advisor:
    def __init__(self, seed: int = 1, rift_at_start: bool = False):
        """`rift_at_start`: the starting plates may rift at once (a Pangaea to break up), instead of
        only after RIFT_TENURE_MYR."""
        self.rng = np.random.default_rng(seed)
        self.rift_at_start = rift_at_start
        self.drift = {}          # plate -> (rotation vector, steps left in its stage)
        self.last_rift = {}      # plate -> time it last rifted (or was made)
        # plate -> time its old ocean was last detached. A clock of its own: on the rift clock, a
        # supercontinent shedding its old ocean every RIDGE_TENURE_MYR never came due to rift.
        self.last_detach = {}
        # A plate made by a rift goes on moving away from its parent for RIFT_TENURE_MYR: child ->
        # (parent, separating rotation, until). Left to drift on its own, it often turned back,
        # closed the young ocean and, with rigid continents, welded again at once.
        self.spreading = {}
        self.pending = {}        # plate proposed to rift -> its separating rotation

    # ------------------------------------------------------------------ proposals
    def propose(self, history) -> tuple:
        """(poles {plate: [lat, lon, rate]}, events [declaration], reasons [str]) for the next step."""
        snapshot = history.current
        world = snapshot.world
        previous = self._previous(history)
        reasons = []
        poles = {}
        pull, perimeter = self._slab_pull(snapshot)
        push = self._ridge_push(world)
        for p in world.plates():
            region = world.regions[p]
            centre = _centre(region)
            continent = sum(c.area_km2() for c in world.crust if c.plate == p and c.kind == 'continent')
            share = continent / max(sphere.area_km2(region), 1.)
            continental = share >= CONTINENTAL_SHARE or continent >= CONTINENTAL_KM2
            torque = pull.get(p, np.zeros(3))
            trench_share = float(np.linalg.norm(torque)) / max(perimeter.get(p, 1.), 1.)
            if np.linalg.norm(torque) > 0.:
                direction = sphere.unit(torque + RIDGE_WEIGHT * push.get(p, np.zeros(3)))
                speed = SLOW_CM_YR + (SLAB_CM_YR - SLOW_CM_YR) * min(3. * trench_share, 1.)
                speed *= 1. - .5 * share
                pulled = True
                self.drift.pop(p, None)
            else:
                direction, speed = self._drifting(p, previous.get(p), push.get(p), continental, centre)
                pulled = False
            if not continental and sphere.area_km2(region) >= LARGE_SHARE * sphere.SPHERE_KM2:
                speed = min(speed, LARGE_CM_YR)   # the global ocean, not a continent's plate
            if pulled:
                reasons.append('plate %d: slab pull, %.1f cm/yr' % (p, min(speed, MAX_CM_YR)))
            rotation = direction * _cm_yr_to_rad_myr(min(speed, MAX_CM_YR), direction, centre)
            old = previous.get(p)
            if old is not None and np.any(old):
                rotation = INERTIA * old + (1. - INERTIA) * rotation
            poles[p] = [round(v, 4) for v in vector_pole(rotation)]
        for child, (parent, away, until) in list(self.spreading.items()):
            if world.time_myr >= until or child not in poles or parent not in poles:
                del self.spreading[child]
                continue
            # It keeps opening its rift, slab or no slab: a piece of a supercontinent is ringed by
            # its parent's trenches from the start, and pulled by them it turned back and welded again.
            rotation = pole_vector(*poles[parent]) + away
            poles[child] = [round(v, 4) for v in vector_pole(rotation)]
        events = self._polarity(snapshot, reasons) + self._rift(world, reasons, poles)
        rifting = {e['plate'] for e in events if e['kind'] == 'rift'}
        events += self._ridge(world, reasons, poles, rifting)
        rifting |= {e['plate'] for e in events if e['kind'] == 'rift'}
        events += self._collapse(world, reasons, poles, rifting)
        return poles, events, reasons

    def _collapse(self, world, reasons, poles, rifting=()) -> list:
        """Detach the old ocean of continental plates that have much of it."""
        events = []
        for p in world.plates():
            if p in rifting or world.time_myr - max(self.last_rift.get(p, 0.), self.last_detach.get(p, 0.)) \
                    < RIDGE_TENURE_MYR:
                continue
            continent = [c for c in world.crust if c.plate == p and c.kind == 'continent']
            if sum(c.area_km2() for c in continent) < CONTINENTAL_KM2:
                continue
            old = [c for c in world.crust if c.plate == p and c.kind == 'ocean' and c.born_myr is not None
                   and world.time_myr - c.born_myr >= OLD_OCEAN_MYR]
            if sum(c.area_km2() for c in old) < OLD_OCEAN_KM2:
                continue
            frame = world.frame(p)
            ocean_at = sphere.rotate(_centre(max(old, key=lambda c: c.area_km2()).shape), frame)
            land_at = sphere.rotate(_centre(max(continent, key=lambda c: c.area_km2()).shape), frame)
            towards = sphere.unit(land_at - np.dot(land_at, ocean_at) * ocean_at)
            axis = sphere.unit(np.cross(ocean_at, towards))
            rate = _cm_yr_to_rad_myr(COLLAPSE_CM_YR, axis, ocean_at)
            parent = pole_vector(*poles[p]) if p in poles else np.zeros(3)
            self.last_detach[p] = world.time_myr
            reasons.append('the old ocean of plate %d is detached, to go down beneath its continent' % p)
            events.append(dict(kind='detach', plate=p, min_age_myr=OLD_OCEAN_MYR,
                               pole=[round(v, 4) for v in vector_pole(parent + axis * rate)]))
        return events

    def _ridge(self, world, reasons, poles, rifting=()) -> list:
        """New mid-ocean ridges across old ocean plates that are due one (at most two a step)."""
        due = []
        for p in world.plates():
            # A plate with a continent is due one too when it has the ocean for one: the plate of a
            # supercontinent with its ocean held 300 Mkm2 of ocean and never had a new ridge.
            area = sum(c.area_km2() for c in world.crust if c.plate == p and c.kind == 'ocean')
            if area < RIDGE_PLATE_SHARE * sphere.SPHERE_KM2:
                continue
            if world.time_myr - self.last_rift.get(p, 0.) < RIDGE_TENURE_MYR:
                continue
            ocean = [(c.area_km2(), world.time_myr - c.born_myr) for c in world.crust
                     if c.plate == p and c.kind == 'ocean' and c.born_myr is not None]
            age = sum(a * t for a, t in ocean) / max(sum(a for a, _ in ocean), 1.)
            if age >= RIDGE_OCEAN_AGE_MYR and p not in rifting:
                due.append((area, p))
        events = []
        for area, plate in sorted(due, reverse=True)[:2]:
            events += self._ridge_one(world, reasons, poles, plate)
        return events

    def _ridge_one(self, world, reasons, poles, plate) -> list:
        # The plate's point farthest from any continent, of a few hundred tried inside it.
        trial = sphere.unit(self.rng.normal(size=(600, 3)))
        lon, lat = sphere.to_lonlat(trial)
        inside = np.asarray(spherely.contains(world.regions[plate], spherely.points(lon, lat)))
        trial = trial[inside]
        land = [sphere.rotate(r, world.frame(c.plate)) for c in world.crust if c.kind == 'continent'
                for poly in sphere.polygons(c.shape) for r in poly[:1]]
        if not len(trial) or not land:
            return []
        land = np.vstack(land)
        far = trial[int(np.argmin(np.max(trial @ land.T, axis=1)))]
        own = [c.shape for c in world.crust if c.plate == plate and c.kind == 'continent']
        own = _main_land(sphere.rotate_geography(sphere.union(*own), world.frame(plate), tolerance=1e-4)) if own else None
        for _ in range(RIFT_TRIES):
            heading = sphere.unit(np.cross(far, sphere.unit(self.rng.normal(size=3))))
            circle = [sphere.unit(far * math.cos(t) + heading * math.sin(t)) for t in np.linspace(-math.pi, math.pi, 73)]
            # A ridge goes through ocean: on a plate with a continent, round a circle clear of it.
            if own is None or _clear_of(own, far, heading) >= RIDGE_CLEAR_KM:
                side = sphere.unit(np.cross(far, heading))             # left of the line, along the sphere
                away = sphere.unit(np.cross(far, side))
                break
        else:
            # No great circle misses the continents (a supercontinent's pieces strewn round its
            # plate): the ocean round the point farthest from them comes away as a plate of its
            # own instead, a ridge opening behind it and a trench ahead.
            far_lon, far_lat = sphere.to_lonlat(far[None])
            clear = float(spherely.distance(spherely.create_point(float(far_lon[0]), float(far_lat[0])), own)) \
                * sphere.RADIUS_KM / spherely.EARTH_RADIUS_METERS - RIDGE_CLEAR_KM
            if clear < RIDGE_CAP_MIN_KM:
                return []
            radius = min(clear, RIDGE_CAP_MAX_KM) / sphere.RADIUS_KM
            east = sphere.unit(np.cross([0., 0., 1.], far) if abs(far[2]) < .99 else np.cross([1., 0., 0.], far))
            north = np.cross(far, east)
            t = np.linspace(0., 2. * math.pi, 73)
            circle = [sphere.unit(math.cos(radius) * far + math.sin(radius) * (math.cos(a) * east + math.sin(a) * north))
                      for a in t]   # anticlockwise seen from outside: the cap is on its left
            away = sphere.unit(np.cross(far, sphere.unit(self.rng.normal(size=3))))   # an axis square to it
        lon, lat = sphere.to_lonlat(np.array(circle))
        rate = _cm_yr_to_rad_myr(4., away, far)
        parent = pole_vector(*poles[plate]) if plate in poles else np.zeros(3)
        self.pending[plate] = sphere.unit(away) * rate
        self.last_rift[plate] = world.time_myr
        reasons.append('a new ridge splits the ocean of plate %d' % plate)
        return [dict(kind='rift', plate=plate, points=[[round(float(a), 3), round(float(b), 3)] for a, b in zip(lon, lat)],
                     new_side='left', pole=[round(v, 4) for v in vector_pole(parent + sphere.unit(away) * rate)])]

    def _shares(self, world) -> dict:
        return {p: sum(c.area_km2() for c in world.crust if c.plate == p and c.kind == 'continent')
                / max(sphere.area_km2(world.regions[p]), 1.) for p in world.plates()}

    def _polarity(self, snapshot, reasons) -> list:
        """Trenches between a continental and an oceanic plate, with the continental one overriding."""
        world = snapshot.world
        share = self._shares(world)
        events = []
        for b in world.boundaries:
            left, right = share.get(b.left, 0.), share.get(b.right, 0.)
            if (left >= CONTINENTAL_SHARE) == (right >= CONTINENTAL_SHARE):
                continue   # both continental, or both oceanic: the core decides
            over = b.left if left > right else b.right
            start = 0
            n = len(b.kinds)
            for k in range(1, n + 1):
                if k < n and (b.kinds[k], b.overriding[k]) == (b.kinds[start], b.overriding[start]):
                    continue
                kind, current = b.kinds[start], b.overriding[start]
                # A collision with an oceanic plate (a terrane brought it) is a trench again.
                wrong = (kind == 'trench' and current != over) or kind in ('converging', 'collision')
                if wrong and k - start >= 3:
                    lon, lat = sphere.to_lonlat(b.points[start:k + 1])
                    events.append(dict(kind='subduction', plates=[b.left, b.right], overriding=over,
                                       points=[[round(float(x), 3), round(float(y), 3)] for x, y in zip(lon, lat)]))
                start = k
        if events:
            reasons.append('%d trench%s with the continental plate overriding' % (len(events), '' if len(events) == 1 else 'es'))
        return events

    def _previous(self, history) -> dict:
        steps = history.project.steps
        if not steps:
            return {}
        last = {int(p): pole_vector(*v) for p, v in steps[-1]['poles'].items()}
        for p, pole in getattr(history.current, 'new_poles', {}).items():
            last[p] = pole_vector(*pole)
        return last

    def _slab_pull(self, snapshot):
        """{plate: torque} from its downgoing trenches, and {plate: perimeter km}."""
        world = snapshot.world
        pull, perimeter = {}, {}
        for b in world.boundaries:
            lengths = sphere.arc_km(b.points[:-1], b.points[1:])
            for p in (b.left, b.right):
                perimeter[p] = perimeter.get(p, 0.) + float(lengths.sum())
            for k, (kind, over) in enumerate(zip(b.kinds, b.overriding)):
                if kind != 'trench' or over not in (b.left, b.right):
                    continue
                down = b.right if over == b.left else b.left
                a, c = b.points[k], b.points[k + 1]
                mid = sphere.unit(a + c)
                left = sphere.unit(np.cross(mid, sphere.unit(c - a)))   # towards b.left, along the sphere
                towards = left if over == b.left else -left
                pull[down] = pull.get(down, np.zeros(3)) + np.cross(mid, towards) * float(lengths[k])
        return pull, perimeter

    def _ridge_push(self, world):
        push = {}
        for b in world.boundaries:
            for k, kind in enumerate(b.kinds):
                if kind != 'ridge':
                    continue
                a, c = b.points[k], b.points[k + 1]
                mid = sphere.unit(a + c)
                length = float(sphere.arc_km(a, c))
                left = sphere.unit(np.cross(mid, sphere.unit(c - a)))
                # Each side is pushed away from the ridge: the left plate to its left.
                push[b.left] = push.get(b.left, np.zeros(3)) + np.cross(mid, left) * length
                push[b.right] = push.get(b.right, np.zeros(3)) + np.cross(mid, -left) * length
        return {p: sphere.unit(v) for p, v in push.items() if np.linalg.norm(v) > 0.}

    def _drifting(self, plate, old, push, continental, centre):
        """A drifting plate's rotation axis and speed. The axis is square to the plate's centre, so the
        plate travels rather than spins: a random axis near the centre turned a continent in place
        at 12 cm/yr round its edges."""
        direction, left = self.drift.get(plate, (None, 0))
        if direction is None or left <= 0:
            base = sphere.unit(old) if old is not None and np.any(old) else sphere.unit(self.rng.normal(size=3))
            turn = sphere.unit(self.rng.normal(size=3))
            direction = base + .6 * turn + (RIDGE_WEIGHT * push if push is not None else 0.)
            direction = sphere.unit(direction - np.dot(direction, centre) * centre)
            # Steer away from the poles: drifting at random, continents came to sit over them.
            lat = math.degrees(math.asin(max(-1., min(1., float(centre[2])))))
            if abs(lat) > POLE_STEER_DEG:
                # About this axis the centre moves towards the equator: (c x z) x c points north.
                equatorward = sphere.unit(np.cross(centre, [0., 0., 1.])) * (-1. if lat > 0 else 1.)
                weight = (abs(lat) - POLE_STEER_DEG) / (90. - POLE_STEER_DEG)
                direction = sphere.unit(direction + 2. * weight * equatorward)
            left = STAGE_STEPS
        self.drift[plate] = (direction, left - 1)
        if continental:
            return direction, self.rng.uniform(*CONTINENTAL_CM_YR)
        return direction, SLOW_CM_YR * self.rng.uniform(.8, 2.)

    def _rift(self, world, reasons, poles) -> list:
        """Rifts of plates carrying a large continent that have not rifted lately, each by chance."""
        candidates = []
        for p in world.plates():
            continent = [c for c in world.crust if c.plate == p and c.kind == 'continent']
            area = sum(c.area_km2() for c in continent)
            # Nothing rifts in its first RIFT_TENURE_MYR, unless the start is a Pangaea to break up.
            since = world.time_myr - self.last_rift.get(p, -np.inf if self.rift_at_start else 0.)
            if area >= RIFT_CONTINENT_KM2 and since >= RIFT_TENURE_MYR:
                candidates.append((area, p, continent))
        events = []
        for area, plate, continent in sorted(candidates, key=lambda c: -c[0]):
            if self.rng.random() <= RIFT_CHANCE:
                events += self._rift_one(world, reasons, poles, area, plate, continent)
        return events

    def _rift_one(self, world, reasons, poles, area, plate, continent) -> list:
        largest = max(continent, key=lambda c: c.area_km2())
        centre = sphere.rotate(_centre(largest.shape), world.frame(plate))
        # Prefer the plate's youngest suture (Scotese XI): through its middle, along it.
        sutures = [s for s in getattr(world, 'sutures', []) if s['plate'] == plate and len(s['points']) >= 2]
        headings = []
        if sutures:
            s = max(sutures, key=lambda s: s['time_myr'])
            points = sphere.rotate(np.asarray(s['points'], float), world.frame(plate))
            centre = sphere.unit(points[len(points) // 2])
            heading = sphere.unit(points[-1] - points[0])
            headings.append(sphere.unit(heading - np.dot(heading, centre) * centre))
        headings += [sphere.unit(np.cross(centre, sphere.unit(self.rng.normal(size=3)))) for _ in range(RIFT_TRIES)]
        land = _main_land(sphere.rotate_geography(sphere.union(*[c.shape for c in continent]), world.frame(plate),
                                                  tolerance=1e-4))
        for heading in headings:
            axis = _opening_axis(land, centre, heading)
            if axis is not None:
                break
        else:
            return []   # every cut tried crosses the continent round more than half its circle
        # A half great circle through the centre, far past the plate's edges (a plate wrapping round
        # its ends is cut along the whole circle: rifting.rift_piece).
        line = [sphere.unit(centre * math.cos(t) + heading * math.sin(t)) for t in np.linspace(-1.5, 1.5, 31)]
        lon, lat = sphere.to_lonlat(np.array(line))
        rate = _cm_yr_to_rad_myr(SEPARATION_CM_YR, axis, centre)
        # The new plate moves with its parent and away from it: a pole of its own in absolute terms
        # let the parent catch it up, collide and weld again within the step.
        parent = pole_vector(*poles[plate]) if plate in poles else np.zeros(3)
        self.pending[plate] = axis * rate
        self.last_rift[plate] = world.time_myr
        reasons.append('plate %d rifts (%.0f million km2 of continent)' % (plate, area / 1e6))
        return [dict(kind='rift', plate=plate, points=[[round(float(a), 3), round(float(b), 3)] for a, b in zip(lon, lat)],
                     new_side='left', pole=[round(v, 4) for v in vector_pole(parent + axis * rate)])]

    def remember(self, history):
        """Recover the rift clocks from a history made before this advisor existed (carrying on a run)."""
        steps = history.project.steps
        for k, snapshot in enumerate(history.snapshots[1:]):
            declared = [e for e in steps[k].get('events', []) if e['kind'] == 'rift'] if k < len(steps) else []
            for e in [e for entry in snapshot.entries for e in entry.get('events', [])]:
                if e.get('kind') == 'plate_split':
                    self._follow(e['new'], e['plate'], snapshot.time_myr)
                if e.get('kind') != 'rift':
                    continue
                self.last_rift[e['new']] = snapshot.time_myr
                if _detached(e):
                    self.last_detach[e['plate']] = snapshot.time_myr
                    continue
                self.last_rift[e['plate']] = snapshot.time_myr
                # The separating rotation is what the rift's pole adds to its parent's.
                rift = next((d for d in declared if d['plate'] == e['plate'] and 'pole' in d), None)
                parent = steps[k]['poles'].get(str(e['plate'])) if rift else None
                if parent is not None:
                    self.spreading[e['new']] = (e['plate'], pole_vector(*rift['pole']) - pole_vector(*parent),
                                                snapshot.time_myr + SEPARATION_MYR)

    def committed(self, snapshot):
        """Note what a committed step made: a plate made by a rift has just rifted, as has its parent."""
        for entry in snapshot.entries[:1]:
            for e in entry.get('events', []):
                if e.get('kind') == 'rift':
                    self.last_rift[e['new']] = snapshot.time_myr
                    if _detached(e):
                        continue   # its old ocean went: not a rift of the plate itself
                    self.last_rift[e['plate']] = snapshot.time_myr
                    if e['plate'] in self.pending:
                        self.spreading[e['new']] = (e['plate'], self.pending[e['plate']],
                                                    snapshot.time_myr + SEPARATION_MYR)
        # A piece a plate came apart into during the step moves on with it, as it did within the
        # step: given a motion of its own, the far piece of a rifted continent turned back and welded.
        for entry in snapshot.entries:
            for e in entry.get('events', []):
                if e.get('kind') == 'plate_split':
                    self._follow(e['new'], e['plate'], snapshot.time_myr)
        self.pending = {}

    def _follow(self, piece, source, time_myr):
        self.last_rift[piece] = time_myr
        if source in self.spreading:
            self.spreading[piece] = self.spreading[source]
        else:
            self.spreading[piece] = (source, np.zeros(3), time_myr + SEPARATION_MYR)
