"""Splitting a plate in two: rifts drawn across a plate, and back-arcs opened behind a trench.

Both make a new plate out of part of an existing one at the start of a step. The new plate starts
in its parent's frame, so its crust keeps its coordinates and it moves with its parent until it is
given a pole of its own. The boundary between them is new, typed from the step's motion: moved
apart, it opens as a ridge.

- A **rift** is a line drawn across a plate, from beyond one edge to beyond another. The side of the
  line the person names (left or right, walking the line as drawn) becomes the new plate.
- A **back-arc** is opened behind a trench. The band of the overriding plate between the trench and
  a line BACKARC_KM inboard of it (the forearc and the arc, ARC_OFFSET_KM in) becomes the new plate.
  Rolled back towards the trench, it leaves a back-arc basin behind it.
"""
from __future__ import annotations

import numpy as np
import spherely

from paradise.core import sphere
from paradise.core.build import densify
from paradise.core.world import Crust

CUT_KM = 1.          # half-width of the cut a rift line makes; its ground stays with the parent
BACKARC_KM = 250.    # inboard of the trench: past the arc (150 km) with room for the basin to open
MINIMUM_PLATE_KM2 = 1e5


def new_plate_id(world) -> int:
    """A plate number never used in this world, the same every time the history is replayed."""
    used = set(world.plates()) | set(world.frames) | set(world.parents)
    return max(used) + 1


def _side(points, probe) -> float:
    """Positive where `probe` lies left of the polyline `points` (seen from outside), at its nearest segment."""
    a, b = points[:-1], points[1:]
    mids = sphere.unit(a + b)
    k = int(np.argmin(sphere.arc_km(mids, probe[None])))
    return float(np.dot(np.cross(a[k], b[k]), probe))


def _centre(geography) -> np.ndarray:
    c = spherely.centroid(geography)
    return sphere.from_lonlat(spherely.get_x(c), spherely.get_y(c))


def rift_piece(world, plate: int, points, new_side: str):
    """The part of `plate` on `new_side` of the line, in the world frame, or ValueError."""
    if plate not in world.regions:
        raise ValueError('there is no plate %d now' % plate)
    line = densify(sphere.from_lonlat(*np.asarray(points, float).T), spacing_km=50.)
    region = world.regions[plate]
    sign = 1. if new_side == 'left' else -1.
    if len(line) > 3 and sphere.arc_km(line[0], line[-1]) < 1.:
        # A closed line (a patch of ocean coming away): the new plate is what it encloses on the
        # named side. Cut as a strip, the strip filled the ring and its inside was lost.
        ring = line[:-1] if sign > 0 else line[:-1][::-1]
        inside = sphere.polygon_robust(ring)
        piece = sphere.intersection(region, inside)
        if sphere.area_km2(piece) < WRAPPED_SHARE * sphere.area_km2(region) or sphere.is_empty(
                sphere.difference(region, piece)):
            raise ValueError('the closed line must enclose part of plate %d, not all of it' % plate)
        return piece
    rest = sphere.difference(region, sphere.strip(line, CUT_KM))
    parts, _ = sphere.parts(rest)
    chosen = [p for p in parts if sign * _side(line, _centre(p)) > 0.]
    kept = [p for p in parts if sign * _side(line, _centre(p)) <= 0.]
    whole = sphere.area_km2(region)
    if not chosen or not kept or min(sphere.area_km2(sphere.union(*chosen)), sphere.area_km2(sphere.union(*kept))) \
            < WRAPPED_SHARE * whole:
        # A long line that does not part the plate goes as far as a line can: the plate wraps round
        # its ends (a supercontinent's plate with its ocean can cover most of the sphere). It is cut
        # along the line's whole great circle. A short one was drawn short, and is refused.
        if float(np.sum(sphere.arc_km(line[:-1], line[1:]))) >= WRAPPED_LINE_KM:
            piece = _hemisphere_piece(region, line, sign)
            if piece is not None:
                return piece
        raise ValueError('the rift line must cross plate %d from one edge to another' % plate)
    return sphere.union(*chosen)


WRAPPED_LINE_KM = 10000.   # a quarter of the way round
WRAPPED_SHARE = .01   # a side with less of the plate than this is a sliver: the line has not parted it


def _hemisphere_piece(region, line, sign):
    """The part of `region` on the `sign` side of the great circle through `line`, if it parts it."""
    normal = sphere.unit(np.sum(np.cross(line[:-1], line[1:]), axis=0)) * sign
    u = sphere.unit(np.cross(normal, line[len(line) // 2]))
    v = np.cross(normal, u)
    t = np.linspace(0., 2. * np.pi, 144, endpoint=False)
    half = sphere.polygon(np.outer(np.cos(t), u) + np.outer(np.sin(t), v))
    lon, lat = sphere.to_lonlat(normal[None])
    if not bool(spherely.contains(half, spherely.create_point(float(lon[0]), float(lat[0])))):
        half = sphere.difference(sphere._full(), half)
    piece = sphere.intersection(region, half)
    area, whole = sphere.area_km2(piece), sphere.area_km2(region)
    if area < WRAPPED_SHARE * whole or area > (1. - WRAPPED_SHARE) * whole:
        return None
    return piece


def backarc_piece(world, trench_points, plates, overriding: int):
    """The band of the overriding plate from the trench to BACKARC_KM inboard, in the world frame."""
    from paradise.terranes import offset_line
    if overriding not in world.regions:
        raise ValueError('there is no plate %d now' % overriding)
    trench = sphere.from_lonlat(*np.asarray(trench_points, float).T)
    trench = trench[np.r_[True, sphere.arc_km(trench[:-1], trench[1:]) > 1e-3]]
    if len(trench) < 2:
        raise ValueError('the trench is too short to open a back-arc behind')
    to_left = overriding == plates[0]
    # Out a little beyond the trench, so the band takes the overriding plate right up to it.
    outer = offset_line(trench, 20., to_left=not to_left)
    inner = offset_line(trench, BACKARC_KM, to_left=to_left)
    ring = np.vstack((outer, inner[::-1]))
    if not to_left:
        ring = ring[::-1]   # the band's interior on the left of its ring
    band = sphere.polygon_robust(ring)
    if sphere.area_km2(band) > .5 * sphere.SPHERE_KM2:
        band = sphere.difference(sphere._full(), band)
    piece = sphere.intersection(band, world.regions[overriding])
    parts, _ = sphere.parts(piece)
    if not parts:
        raise ValueError('there is no room for a back-arc behind this trench')
    return max(parts, key=sphere.area_km2)


EDGE_KM = 1.   # a boundary segment this close to a split-off piece runs along its edge


def _split_by_piece(boundaries, plate, new, piece):
    """The boundaries with `plate`'s segments along the edge of `piece` given to `new`."""
    from paradise.core.regions import Boundary, _relabel
    result = []
    for b in boundaries:
        if plate not in (b.left, b.right) or len(b.points) < 2:
            result.append(b)
            continue
        lon, lat = sphere.to_lonlat(b.mids())
        km = np.asarray(spherely.distance(spherely.points(lon, lat), piece)) * sphere.RADIUS_KM / spherely.EARTH_RADIUS_METERS
        side = km <= EDGE_KM
        start = 0
        for k in range(1, len(side) + 1):
            if k < len(side) and side[k] == side[start]:
                continue
            idents = None if getattr(b, 'idents', None) is None else np.asarray(b.idents)[start:k]
            part = Boundary(b.left, b.right, b.points[start:k + 1], b.kinds[start:k], b.overriding[start:k],
                            b.closed_km[start:k].copy(), idents)
            result += _relabel([part], {plate: new} if side[start] else {})
            start = k
    return result


def old_ocean_piece(world, plate: int, min_age_myr: float):
    """The largest connected stretch of `plate`'s ocean older than `min_age_myr`, in the world frame."""
    if plate not in world.regions:
        raise ValueError('there is no plate %d now' % plate)
    old = [c.shape for c in world.crust if c.plate == plate and c.kind == 'ocean' and c.born_myr is not None
           and world.time_myr - c.born_myr >= min_age_myr]
    if not old:
        raise ValueError('plate %d has no ocean %.0f Myr old' % (plate, min_age_myr))
    shape = sphere.rotate_geography(sphere.union(*old), world.frame(plate), tolerance=1e-4)
    parts, _ = sphere.parts(sphere.intersection(shape, world.regions[plate]))
    if not parts:
        raise ValueError('plate %d has no ocean %.0f Myr old' % (plate, min_age_myr))
    return max(parts, key=sphere.area_km2)


def split_off(world, plate: int, piece, rotations: dict, name: str, cratons: list, features: list,
              rifts: list, line=None) -> tuple:
    """Make `piece` (world frame) of `plate` a new plate; returns (new plate, cratons, event).

    The crust, the cratons, the arcs and terranes and the sutures on the piece go with it.
    """
    area = sphere.area_km2(piece)
    whole = sphere.area_km2(world.regions[plate])
    if area < MINIMUM_PLATE_KM2 or whole - area < MINIMUM_PLATE_KM2:
        raise ValueError('splitting plate %d there leaves a piece smaller than %.0f thousand km2'
                         % (plate, MINIMUM_PLATE_KM2 / 1e3))
    new = new_plate_id(world)
    frame = world.frame(plate)
    world.frames[new] = frame.copy()
    world.parents[new] = plate
    world.names[new] = name
    local = sphere.rotate_geography(piece, frame.T, tolerance=1e-4)
    moved = []
    for c in world.crust:
        if c.plate != plate:
            continue
        taken = sphere.intersection(c.shape, local)
        if sphere.area_km2(taken) >= sphere.DETAIL_KM2:
            c.shape = sphere.difference(c.shape, taken)
            # Same frame, so the same coordinates; a continent's thickness mesh is rebuilt from its mean.
            moved.append(Crust(new, c.kind, taken, c.born_myr, c.thickness_km))
    world.crust = [c for c in world.crust if sphere.area_km2(c.shape) >= sphere.DETAIL_KM2] + moved

    def on_piece(shape):
        return spherely.contains(local, spherely.centroid(shape))

    cratons = [(new, shape) if p == plate and on_piece(shape) else (p, shape) for p, shape in cratons]
    for f in features:
        if f.plate == plate and on_piece(f.shape):
            f.plate = new
    for s in getattr(world, 'sutures', []):
        if s['plate'] == plate:
            lon, lat = sphere.to_lonlat(sphere.unit(np.mean(s['points'], axis=0)))
            if spherely.contains(local, spherely.create_point(float(lon), float(lat))):
                s['plate'] = new
    for r in rifts:
        if r['plate'] == plate and len(r['points']) and spherely.contains(
                local, spherely.create_point(*map(float, sphere.to_lonlat(sphere.unit(np.mean(r['points'], axis=0)))))):
            r['plate'] = new
    # The parent's boundaries along the piece's edge now belong to the new plate: a trench the
    # forearc takes over stays a trench, with the forearc overriding it. The core only carries a
    # boundary's state between the same two plates.
    old = _split_by_piece(world.boundaries, plate, new, piece)
    world.regions = world._from_crust()
    full = {p: rotations.get(p, np.zeros(3)) for p in world.plates()}
    world.boundaries = world._derive(world._index(old), full)
    if line is not None:
        rifts.append(dict(plate=plate, new=new, time_myr=world.time_myr, kind=name,
                          points=sphere.rotate(line, frame.T)))
    return new, cratons, dict(kind='rift', plate=plate, new=new, area_km2=round(area), name=name)
