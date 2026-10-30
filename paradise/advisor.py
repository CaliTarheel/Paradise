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
RIDGE_WEIGHT = .25
INERTIA = .5             # how much of a plate's old rotation each new one keeps
STAGE_STEPS = 5
RIFT_CONTINENT_KM2 = 10e6   # Artifexia's continents are mostly 10 to 50 million km2
RIFT_TENURE_MYR = 60.      # Artifexia's continents split every 50 to 100 Myr
RIFT_CHANCE = .35        # per step, once a plate may rift
# A large ocean plate splits along a new mid-ocean ridge, as a Pasta world adds ridges to its one
# ocean plate: without it, ocean leaves only where continents run over it, and grows ever older.
# Artifexia has 30 to 46 plates a stage: ocean plates are split far more often than one per globe.
RIDGE_PLATE_SHARE, RIDGE_OCEAN_AGE_MYR, RIDGE_TENURE_MYR = .08, 80., 40.


def _cm_yr_to_rad_myr(cm_yr, axis, centre):
    """The rotation rate (radians per Myr) about `axis` that moves `centre` at `cm_yr`."""
    s = float(np.linalg.norm(np.cross(axis, centre)))
    return (cm_yr * 10. / sphere.RADIUS_KM) / max(s, .2)


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
            torque = pull.get(p, np.zeros(3))
            trench_share = float(np.linalg.norm(torque)) / max(perimeter.get(p, 1.), 1.)
            if np.linalg.norm(torque) > 0.:
                direction = sphere.unit(torque + RIDGE_WEIGHT * push.get(p, np.zeros(3)))
                speed = SLOW_CM_YR + (SLAB_CM_YR - SLOW_CM_YR) * min(3. * trench_share, 1.)
                speed *= 1. - .5 * share
                pulled = True
                self.drift.pop(p, None)
            else:
                direction, speed = self._drifting(p, previous.get(p), push.get(p), share >= CONTINENTAL_SHARE)
                pulled = False
            if sphere.area_km2(region) >= LARGE_SHARE * sphere.SPHERE_KM2:
                speed = min(speed, LARGE_CM_YR)
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
            if child not in pull:   # unless a slab now pulls it, it keeps opening its rift
                rotation = pole_vector(*poles[parent]) + away
                poles[child] = [round(v, 4) for v in vector_pole(rotation)]
        events = self._polarity(snapshot, reasons) + self._rift(world, reasons, poles)
        rifting = {e['plate'] for e in events if e['kind'] == 'rift'}
        events += self._ridge(world, reasons, poles, rifting)
        return poles, events, reasons

    def _ridge(self, world, reasons, poles, rifting=()) -> list:
        """New mid-ocean ridges across old ocean plates that are due one (at most two a step)."""
        shares = self._shares(world)
        due = []
        for p in world.plates():
            area = sphere.area_km2(world.regions[p])
            if shares[p] >= .1 or area < RIDGE_PLATE_SHARE * sphere.SPHERE_KM2:
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
        heading = sphere.unit(np.cross(far, sphere.unit(self.rng.normal(size=3))))
        circle = [sphere.unit(far * math.cos(t) + heading * math.sin(t)) for t in np.linspace(-math.pi, math.pi, 73)]
        lon, lat = sphere.to_lonlat(np.array(circle))
        side = sphere.unit(np.cross(far, heading))             # left of the line, along the sphere
        away = np.cross(far, side)
        rate = _cm_yr_to_rad_myr(4., sphere.unit(away), far)
        parent = pole_vector(*poles[plate]) if plate in poles else np.zeros(3)
        self.pending[plate] = sphere.unit(away) * rate
        self.last_rift[plate] = world.time_myr
        reasons.append('a new ridge splits ocean plate %d' % plate)
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
                wrong = (kind == 'trench' and current != over) or kind == 'converging'
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

    def _drifting(self, plate, old, push, continental=False):
        direction, left = self.drift.get(plate, (None, 0))
        if direction is None or left <= 0:
            base = sphere.unit(old) if old is not None and np.any(old) else sphere.unit(self.rng.normal(size=3))
            turn = sphere.unit(self.rng.normal(size=3))
            direction = sphere.unit(base + .6 * turn + (RIDGE_WEIGHT * push if push is not None else 0.))
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
        if sutures:
            s = max(sutures, key=lambda s: s['time_myr'])
            points = sphere.rotate(np.asarray(s['points'], float), world.frame(plate))
            centre = sphere.unit(points[len(points) // 2])
            heading = sphere.unit(points[-1] - points[0])
            heading = sphere.unit(heading - np.dot(heading, centre) * centre)
        else:
            heading = sphere.unit(np.cross(centre, sphere.unit(self.rng.normal(size=3))))
        # A half great circle through the centre, far past the plate's edges.
        line = [sphere.unit(centre * math.cos(t) + heading * math.sin(t)) for t in np.linspace(-1.5, 1.5, 31)]
        lon, lat = sphere.to_lonlat(np.array(line))
        side = sphere.unit(np.cross(centre, heading))          # left of the line, along the sphere
        away = np.cross(centre, side)                           # turns the left side away from the line
        rate = _cm_yr_to_rad_myr(3., sphere.unit(away), centre)
        # The new plate moves with its parent and away from it: a pole of its own in absolute terms
        # let the parent catch it up, collide and weld again within the step.
        parent = pole_vector(*poles[plate]) if plate in poles else np.zeros(3)
        self.pending[plate] = sphere.unit(away) * rate
        self.last_rift[plate] = world.time_myr
        reasons.append('plate %d rifts (%.0f million km2 of continent)' % (plate, area / 1e6))
        return [dict(kind='rift', plate=plate, points=[[round(float(a), 3), round(float(b), 3)] for a, b in zip(lon, lat)],
                     new_side='left', pole=[round(v, 4) for v in vector_pole(parent + sphere.unit(away) * rate)])]

    def committed(self, snapshot):
        """Note what a committed step made: a plate made by a rift has just rifted, as has its parent."""
        for entry in snapshot.entries[:1]:
            for e in entry.get('events', []):
                if e.get('kind') == 'rift':
                    self.last_rift[e['new']] = snapshot.time_myr
                    self.last_rift[e['plate']] = snapshot.time_myr
                    if e['plate'] in self.pending:
                        self.spreading[e['new']] = (e['plate'], self.pending[e['plate']],
                                                    snapshot.time_myr + RIFT_TENURE_MYR)
        self.pending = {}
