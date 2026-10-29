"""A painted starting world: equirectangular grids of plate and crust, turned into plates on the sphere.

The editor paints two grids of the same size, rows from the north pole down and columns from
180 W eastwards:

- `plates`: the plate number of each cell (1 or more);
- `crust`: 0 for ocean, 1 for continent, 2 for craton (continent that is also a craton).

Each (plate, ocean or continent) piece is traced along cell edges into rings whose interior lies to
the left, which is what the core's oriented polygons need, so a plate larger than a hemisphere is
still the region it encloses. The rings of every piece are simplified together
(`sphere.simplify_shared`), so neighbouring pieces still share their edges and the plates still tile
the sphere, with the staircase of the cells taken out.
"""
from __future__ import annotations

import math
import zlib
from collections import defaultdict

import numpy as np
from scipy import ndimage

from paradise.core import sphere
from paradise.core.regions import FULL, RegionWorld
from paradise.core.world import Crust

OCEAN, CONTINENT, CRATON = 0, 1, 2
# Coastlines (edges between a plate's continent and its ocean) get the irregular shape of a
# hand-drawn coast instead of the painting's smoothed cells: a smooth displacement of up to about
# DETAIL_KM on features about a thousand km across, then midpoint displacement down to
# DETAIL_SPACING_KM. The target is the hand-built Artifexia project's level of detail, whose
# continents have a vertex every 270 km (median; 100 to 700 km), not a finer one: shape, not
# vertex count, was what the painted coasts lacked. Plate boundaries take the smooth displacement
# too, so they lose the painting's ruler-straight edges, but no midpoint detail.
# 300 km is about as far as the smooth displacement goes before rings cross themselves (450 km
# crossed on a nine-plate world, and then the coast keeps its painted shape).
DETAIL_KM, DETAIL_SPACING_KM = 300., 100.
# Midpoint displacement, as a fraction of the edge's length either way. Much above this a sharp
# turn's children cross near its tip, and the ring keeps its plain edges instead.
ROUGHNESS = .25
# The longest edge simplification may leave, in cells. A declaration (a trench, say) applies to
# whole boundary segments, so this is also how closely a declared stretch can be placed.
EDGE_CELLS = 2.


def cell_km(width: int) -> float:
    """Width of a cell at the equator."""
    return 2. * math.pi * sphere.RADIUS_KM / width


def corner_points(height: int, width: int) -> np.ndarray:
    """Unit vectors of the cell corners, (height + 1, width): row 0 is the north pole."""
    lat = 90. - 180. * np.arange(height + 1) / height
    lon = -180. + 360. * np.arange(width) / width
    lon, lat = np.meshgrid(lon, lat)
    return sphere.from_lonlat(lon, lat)


def _components(mask: np.ndarray) -> tuple:
    """4-connected components of `mask`, joined across the 180 degree meridian."""
    labels, count = ndimage.label(mask)
    if count == 0:
        return labels, 0
    # Join components that meet across the date line.
    parent = list(range(count + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in zip(labels[:, 0], labels[:, -1]):
        if a and b:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    roots = sorted({find(a) for a in range(1, count + 1)})
    renumber = np.zeros(count + 1, int)
    for a in range(1, count + 1):
        renumber[a] = roots.index(find(a)) + 1
    return renumber[labels], len(roots)


def trace_rings(mask: np.ndarray) -> list:
    """Rings round the cells of `mask` with the cells on their left, as lists of (row, column) corners.

    Columns wrap round the sphere. Edges along the poles are kept (both of their corners are the
    pole); they leave repeated points that `_ring_points` takes out. Where two cells of the mask meet
    only at a corner, the ring turns left, round the cell it came along, so they stay apart.
    """
    height, width = mask.shape
    up = np.roll(mask, 1, axis=0)
    up[0] = False
    down = np.roll(mask, -1, axis=0)
    down[-1] = False
    left = np.roll(mask, 1, axis=1)
    right = np.roll(mask, -1, axis=1)
    edges = defaultdict(list)   # start corner -> [end corner]
    rows, cols = np.nonzero(mask)
    for i, j in zip(rows.tolist(), cols.tolist()):
        j1 = (j + 1) % width
        # Counter-clockwise round the cell seen from outside (east to the right, north up).
        if not down[i, j]:
            edges[(i + 1, j)].append((i + 1, j1))
        if not right[i, j]:
            edges[(i + 1, j1)].append((i, j1))
        if not up[i, j]:
            edges[(i, j1)].append((i, j))
        if not left[i, j]:
            edges[(i, j)].append((i + 1, j))

    def direction(a, b):
        di = b[0] - a[0]
        dj = (b[1] - a[1] + width // 2) % width - width // 2
        return di, dj

    rings = []
    while edges:
        start = next(iter(edges))
        ring, here, came = [], start, None
        while True:
            ring.append(here)
            options = edges[here]
            k = 0
            if len(options) > 1 and came is not None:
                # Turn left: in (row, column) terms, left of a heading (di, dj) is (-dj, di).
                want = (-came[1], came[0])
                k = next((k for k, o in enumerate(options) if direction(here, o) == want), 0)
            nxt = options.pop(k)
            if not options:
                del edges[here]
            came = direction(here, nxt)
            here = nxt
            if here == start:
                break
        rings.append(ring)
    return rings


def _ring_points(ring, corners) -> np.ndarray | None:
    """A traced ring as unit vectors, with repeated points (at the poles) taken out."""
    points = np.array([corners[i, j] for i, j in ring])
    keep = np.ones(len(points), bool)
    for k in range(1, len(points)):
        keep[k] = not np.allclose(points[k], points[k - 1], atol=1e-12)
    points = points[keep]
    if len(points) > 1 and np.allclose(points[0], points[-1], atol=1e-12):
        points = points[:-1]
    if len(points) < 3:
        return None
    return points


def _valid(ring) -> bool:
    try:
        sphere.polygon(ring)
        return True
    except (ValueError, RuntimeError):
        return False


def _key(point):
    return tuple(np.round(point, 9))


def _edge_key(a, b):
    ka, kb = _key(a), _key(b)
    return (ka, kb) if ka < kb else (kb, ka)


def displacement(points, km, seed=1) -> np.ndarray:
    """A fractal displacement field: each point moved up to about `km` along the sphere.

    It depends on position alone, so a vertex two rings share moves the same in both, and the
    pieces still tile the sphere.
    """
    from paradise.relief import fbm
    points = np.asarray(points, float)
    field = np.stack([fbm(points, octaves=4, frequency=3., seed=seed + 31 * k, persistence=.55) for k in range(3)],
                     axis=-1)
    field -= np.sum(field * points, axis=-1, keepdims=True) * points   # along the sphere only
    return sphere.unit(points + field * (2. * km / sphere.RADIUS_KM))


def _midpoints(first, last, seed) -> list:
    """Points between `first` and `last` by recursive midpoint displacement, from `first`.

    Each midpoint moves square to its edge by up to ROUGHNESS of the edge's length, with a
    random number drawn from the two ends alone, so the two rings sharing an edge get the very
    same points; halving goes on until edges are DETAIL_SPACING_KM or less.
    """
    length = float(sphere.arc_km(first, last))
    if length <= DETAIL_SPACING_KM:
        return []
    middle = sphere.unit(first + last)
    side = sphere.unit(np.cross(first, last))
    draw = zlib.crc32(np.round(np.concatenate((first, last)), 9).tobytes() + seed.to_bytes(4, 'little'))
    r = draw / 2. ** 31 - 1.
    middle = sphere.unit(middle + side * (r * ROUGHNESS * length / sphere.RADIUS_KM))
    return _midpoints(first, middle, seed) + [middle] + _midpoints(middle, last, seed)


def _add_detail(ring, coast, km, seed):
    """`ring` with fractal detail on its coast edges.

    Every vertex is moved by the smooth `displacement` field (up to about `km`), then each coast
    edge is divided by midpoint displacement. An edge is always divided from its smaller end,
    so the two rings sharing it get the very same points and the pieces still tile the sphere.
    """
    n = len(ring)
    on_coast = np.array([_edge_key(ring[k], ring[(k + 1) % n]) in coast for k in range(n)])
    # Every vertex moves, plate boundaries' too: moving only the coast's would bend the boundary
    # edges that meet a coast sharply enough to cross. The field depends on position alone, so
    # shared vertices stay shared.
    points = displacement(np.array(ring, float), km, seed)
    out = []
    for k in range(n):
        a, b = points[k], points[(k + 1) % n]
        out.append(a)
        if on_coast[k]:
            forward = _key(a) < _key(b)
            inner = _midpoints(a, b, seed) if forward else _midpoints(b, a, seed)[::-1]
            out += inner
    return np.array(out)


def pieces_from_grid(plates: np.ndarray, kinds: np.ndarray, simplify_km2: float | None = None,
                     detail_km: float = 0., seed: int = 1) -> list:
    """[(plate, kind, geography)] for every connected piece of each (plate, kind) in the grids.

    `kinds` holds 'ocean' or 'continent' as 0 and 1. Each piece is the region left of all of its rings
    (on a sphere every ring of a connected piece has the piece on the same side).
    """
    height, width = plates.shape
    corners = corner_points(height, width)
    key = plates.astype(np.int64) * 2 + kinds.astype(np.int64)
    traced = []   # (plate, kind, component, ring points)
    for value in np.unique(key):
        plate, kind = int(value // 2), int(value % 2)
        labels, count = _components(key == value)
        for c in range(1, count + 1):
            rings = [_ring_points(r, corners) for r in trace_rings(labels == c)]
            rings = [r for r in rings if r is not None]
            traced.append((plate, kind, c, rings))
    if simplify_km2 is None:
        simplify_km2 = 1.2 * cell_km(width) ** 2
    flat, owner = [], []
    for n, (_, _, _, rings) in enumerate(traced):
        for r in rings:
            flat.append(r)
            owner.append(n)
    if flat:
        flat, _, _ = sphere.simplify_shared(flat, area_km2=simplify_km2, max_edge_km=EDGE_CELLS * cell_km(width),
                                            groups=owner)
    if detail_km and flat:
        coast = {_edge_key(a, b) for n, ring in zip(owner, flat) if traced[n][1] == 1
                 for a, b in zip(ring, np.roll(ring, -1, axis=0))}
        plain = flat
        for _ in range(8):
            detailed = [_add_detail(ring, coast, detail_km, seed) for ring in plain]
            # Where detail makes a ring cross itself (a sharp corner of the painting, say), that
            # ring's coast keeps its plain edges, in every ring that shares them, and the rest is
            # tried again.
            crossed = [k for k, ring in enumerate(detailed) if not _valid(ring)]
            if not crossed:
                break
            for k in crossed:
                ring = plain[k]
                coast -= {_edge_key(a, b) for a, b in zip(ring, np.roll(ring, -1, axis=0))}
        else:
            detailed = plain
        flat = detailed
    by_piece = defaultdict(list)
    for n, ring in zip(owner, flat):
        if len(ring) >= 3:
            by_piece[n].append(ring)
    pieces = []
    for n, (plate, kind, _, _) in enumerate(traced):
        rings = by_piece.get(n, [])
        if not rings:
            region = FULL   # one piece covering the whole sphere
        else:
            region = sphere.polygon(rings[0])
            for ring in rings[1:]:
                region = sphere.intersection(region, sphere.polygon(ring))
        if not sphere.is_empty(region):
            pieces.append((plate, 'continent' if kind else 'ocean', region))
    return pieces


def world_from_grid(plates: np.ndarray, crust: np.ndarray, ocean_age_myr: float = 60.,
                    names: dict | None = None, detail_km: float = 0., seed: int = 1) -> tuple:
    """A RegionWorld at time 0 from painted grids, and the cratons as [(plate, geography)].

    The ocean starts `ocean_age_myr` old everywhere. Cratons are continent in the core (they move and
    deform with it); their outlines are kept beside the world, in their plate's frame, so they can
    be exported and later made harder to rift.
    """
    plates = np.asarray(plates, np.int64)
    crust = np.asarray(crust, np.int64)
    if plates.shape != crust.shape:
        raise ValueError('the plate and crust grids differ in size')
    if plates.min() < 1:
        raise ValueError('every cell needs a plate (numbers from 1)')
    pieces = pieces_from_grid(plates, (crust > 0).astype(np.int64), detail_km=detail_km, seed=seed)
    world_crust = [Crust(p, kind, shape, born_myr=-ocean_age_myr if kind == 'ocean' else None)
                   for p, kind, shape in pieces]
    world = RegionWorld(world_crust, {}, {}, [], names=dict(names or {}))
    for c in world.crust:
        world.frame(c.plate)
    world.regions = world._from_crust()
    world.boundaries = world._derive({}, {})
    cratons = []
    if (crust == CRATON).any():
        craton_plates = np.where(crust == CRATON, plates, 0)
        for value in np.unique(craton_plates[craton_plates > 0]):
            mask = craton_plates == value
            labels, count = _components(mask)
            corners = corner_points(*plates.shape)
            for c in range(1, count + 1):
                rings = [_ring_points(r, corners) for r in trace_rings(labels == c)]
                rings = [r for r in rings if r is not None]
                if not rings:
                    continue
                rings, _, _ = sphere.simplify_shared(rings, area_km2=1.2 * cell_km(plates.shape[1]) ** 2,
                                                     max_edge_km=EDGE_CELLS * cell_km(plates.shape[1]))
                region = sphere.polygon(rings[0])
                for ring in rings[1:]:
                    region = sphere.intersection(region, sphere.polygon(ring))
                cratons.append((int(value), region))
    return world, cratons
