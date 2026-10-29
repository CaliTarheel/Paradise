"""Mountain belts, filed by age as a hand-built reconstruction files them: active, former and old.

Belts are derived from the world at the end of each step, not carried:

- **active**: continental crust within HALF_WIDTH_KM of a boundary where continents collide now;
- **former**: continental crust within HALF_WIDTH_KM of a suture made less than FORMER_MYR ago,
  where two plates welded after a collision, and the mountains are still high;
- **old**: the same along older sutures, worn down but still a belt of deformed crust (and a
  weakness a later rift may follow: Scotese's sutures are long-lived).

Each belt is kept in its plate's frame, like every other feature, so the export can carry it.
"""
from __future__ import annotations

import numpy as np

from paradise.core import sphere

HALF_WIDTH_KM = 250.
FORMER_MYR = 100.
MINIMUM_KM2 = 1000.


def _continent(world, plate):
    shapes = [c.shape for c in world.crust if c.plate == plate and c.kind == 'continent']
    return sphere.union(*shapes) if shapes else None


def _belt(world, plate, line_world, kind, since_myr):
    continent = _continent(world, plate)
    if continent is None or len(line_world) < 2:
        return None
    zone = sphere.rotate_geography(sphere.strip(line_world, HALF_WIDTH_KM), world.frame(plate).T, tolerance=1e-4)
    shape = sphere.intersection(continent, zone)
    if sphere.area_km2(shape) < MINIMUM_KM2:
        return None
    return dict(kind=kind, plate=plate, shape=shape, since_myr=since_myr)


def belts(world) -> list:
    """[{kind, plate, shape (plate frame), since_myr}] for every belt in `world` now."""
    result = []
    for b in world.boundaries:
        n = len(b.kinds)
        start = 0
        for k in range(1, n + 1):
            if k < n and b.kinds[k] == b.kinds[start]:
                continue
            if b.kinds[start] == 'collision':
                line = b.points[start:k + 1]
                for plate in (b.left, b.right):
                    belt = _belt(world, plate, line, 'active', world.time_myr)
                    if belt:
                        result.append(belt)
            start = k
    for s in getattr(world, 'sutures', []):
        if s['plate'] not in world.regions:
            continue
        age = world.time_myr - s['time_myr']
        line = sphere.rotate(np.asarray(s['points'], float), world.frame(s['plate']))
        belt = _belt(world, s['plate'], line, 'former' if age < FORMER_MYR else 'old', s['time_myr'])
        if belt:
            result.append(belt)
    return result
