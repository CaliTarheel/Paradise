"""Starting worlds made by rule instead of by hand: plates, a supercontinent, cratons, microcontinents.

The shapes are what an author might sketch, not noise. Plates are a weighted (power) Voronoi
diagram of the sphere, so their sizes vary as Earth's do. Each continental block is the union of
a few overlapping ellipses, which gives the bays, capes and lobes of a drawn outline without any
fractal roughness. The blocks are packed round one centre into a supercontinent, each with a craton
at its core, and a few microcontinents lie offshore.
"""
from __future__ import annotations

import math

import numpy as np

from paradise.core import sphere
from paradise.paint import CONTINENT, CRATON

H, W = 360, 720


def grid_xyz():
    lat = 90. - (np.arange(H) + .5) * 180. / H
    lon = -180. + (np.arange(W) + .5) * 360. / W
    lon, lat = np.meshgrid(lon, lat)
    return sphere.from_lonlat(lon, lat)


def _frame(centre):
    """East and north unit vectors at `centre`."""
    east = np.cross([0., 0., 1.], centre)
    if np.linalg.norm(east) < 1e-9:
        east = np.array([0., 1., 0.])
    east = sphere.unit(east)
    return east, np.cross(centre, east)


def ellipse(xyz, centre, a_km, b_km, angle_deg):
    """Cells inside an ellipse of semi-axes `a_km` and `b_km`, turned `angle_deg` from east."""
    east, north = _frame(centre)
    t = math.radians(angle_deg)
    u = math.cos(t) * east + math.sin(t) * north
    v = -math.sin(t) * east + math.cos(t) * north
    x = np.arctan2(xyz @ u, xyz @ centre) * sphere.RADIUS_KM
    y = np.arctan2(xyz @ v, xyz @ centre) * sphere.RADIUS_KM
    return (xyz @ centre > 0.) & ((x / a_km) ** 2 + (y / b_km) ** 2 <= 1.)


def offset(centre, km, bearing_deg):
    """The point `km` from `centre` along `bearing_deg` (clockwise from north)."""
    east, north = _frame(centre)
    b = math.radians(bearing_deg)
    direction = math.sin(b) * east + math.cos(b) * north
    angle = km / sphere.RADIUS_KM
    return sphere.unit(centre * math.cos(angle) + direction * math.sin(angle))


def block(xyz, centre, size_km, rng):
    """A continental block: a main ellipse and three or four lobes, and a craton in its middle."""
    land = ellipse(xyz, centre, size_km, size_km * rng.uniform(.55, .8), rng.uniform(0., 180.))
    for _ in range(int(rng.integers(3, 5))):
        lobe = offset(centre, size_km * rng.uniform(.5, .9), rng.uniform(0., 360.))
        land |= ellipse(xyz, lobe, size_km * rng.uniform(.35, .6), size_km * rng.uniform(.2, .4),
                        rng.uniform(0., 180.))
    craton = ellipse(xyz, centre, size_km * .4, size_km * .3, rng.uniform(0., 180.))
    return land, craton


def supercontinent_world(seed=7, plates=30, blocks=8, microcontinents=5, centre=(20., 0.)):
    """(plate grid, crust grid) for a supercontinent world."""
    rng = np.random.default_rng(seed)
    xyz = grid_xyz()
    # Plates: power Voronoi, a weight per seed, so some plates are large and some small.
    seeds = sphere.unit(rng.normal(size=(plates, 3)))
    weights = rng.uniform(0., .12, size=plates)
    plate = 1 + np.argmax(xyz @ seeds.T + weights, axis=-1)
    crust = np.zeros(plate.shape, np.int64)
    middle = sphere.from_lonlat(*centre)
    # Blocks round the centre: one in the middle, the rest in a ring, touching or nearly.
    places = [middle] + [offset(middle, rng.uniform(4300., 5000.), 360. * k / (blocks - 1) + rng.uniform(-12., 12.))
                         for k in range(blocks - 1)]
    for place in places:
        land, craton = block(xyz, place, rng.uniform(2300., 2900.), rng)
        crust[land & (crust == 0)] = CONTINENT
        crust[craton] = CRATON
    # Microcontinents well offshore.
    for _ in range(microcontinents):
        place = offset(middle, rng.uniform(8500., 11000.), rng.uniform(0., 360.))
        land = ellipse(xyz, place, rng.uniform(300., 600.), rng.uniform(150., 300.), rng.uniform(0., 180.))
        crust[land & (crust == 0)] = CONTINENT
    return plate, crust


def pasta_world(seed=3, continents=3, margin_km=500., separation_km=8500.):
    """(plate grid, crust grid) as Worldbuilding Pasta starts a world: one ocean plate round the globe
    and a plate per continent, each its continent and a passive margin of ocean `margin_km` wide.

    Nothing moves until the continental plates are given poles: subduction starts where they push
    into the ocean plate, and ridges open behind them.
    """
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    xyz = grid_xyz().reshape(-1, 3)
    centres = []
    while len(centres) < continents:
        c = sphere.unit(rng.normal(size=3))
        c = sphere.unit(c * [1., 1., .6])   # keep away from the poles
        if all(float(sphere.arc_km(c, o)) >= separation_km for o in centres):
            centres.append(c)
    crust = np.zeros(len(xyz), np.int64)
    owner = np.zeros(len(xyz), np.int64)
    for k, c in enumerate(centres):
        land, craton = block(xyz, c, rng.uniform(3700., 4200.), rng)
        land &= owner == 0   # first come, first served where two would touch
        crust[land & (crust == 0)] = CONTINENT
        crust[craton & land] = CRATON
        owner[land] = k + 2
    plate = np.ones(len(xyz), np.int64)
    land = np.flatnonzero(owner)
    tree = cKDTree(xyz[land])
    chord, nearest = tree.query(xyz)
    within = 2. * np.arcsin(np.clip(chord / 2., 0., 1.)) * sphere.RADIUS_KM <= margin_km
    plate[within] = owner[land][nearest[within]]
    return plate.reshape(H, W), crust.reshape(H, W)


def pangaea_world(seed=7, blocks=8, margin_km=500., centre=(20., 0.)):
    """(plate grid, crust grid) as the hand-built Artifexia project starts: one supercontinent on one
    plate, in one ocean plate round the globe.

    The supercontinent is `blocks` continental blocks packed round a centre, each with a craton; its
    plate is the continent and a passive margin of ocean `margin_km` wide. It is all there is to
    break up: rifting it is the first thing that happens.
    """
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    xyz = grid_xyz()
    crust = np.zeros(xyz.shape[:2], np.int64)
    middle = sphere.from_lonlat(*centre)
    places = [middle] + [offset(middle, rng.uniform(4300., 5000.), 360. * k / (blocks - 1) + rng.uniform(-12., 12.))
                         for k in range(blocks - 1)]
    for place in places:
        land, craton = block(xyz, place, rng.uniform(2550., 3100.), rng)
        crust[land & (crust == 0)] = CONTINENT
        crust[craton] = CRATON
    flat = xyz.reshape(-1, 3)
    land = np.flatnonzero(crust.ravel() > 0)
    chord, _ = cKDTree(flat[land]).query(flat)
    within = 2. * np.arcsin(np.clip(chord / 2., 0., 1.)) * sphere.RADIUS_KM <= margin_km
    plate = np.where(within, 2, 1).reshape(crust.shape)
    return plate, crust
