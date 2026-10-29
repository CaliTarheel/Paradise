"""Mantle plumes: hotspots that leave trails, and large igneous provinces (LIPs) where a plume head arrives.

Both are declared, at a point on the map at the start of a step.

- A **hotspot** stays at its point, fixed in the mantle (the anchor plate's frame). After every
  sub-step the plate over it gets a trail mark, kept in that plate's frame, so the trail rides
  away with the plate. Scotese: hotspots give a reference frame, and do not break continents.
- A **LIP** is a plume head: a cap of PLUME_RADIUS_KM (or as declared) emplaced at once.
  - Over continent it is a flood-basalt province, recorded as a feature.
  - Over ocean the crust becomes an oceanic plateau: thick, buoyant crust, made continental in the
    core so that it resists subduction and can later accrete as a terrane (Ontong Java, say).
  A LIP leaves a hotspot at its centre, the plume's tail, unless declared without one.
"""
from __future__ import annotations

import math

import numpy as np
import spherely

from paradise.core import sphere
from paradise.core.world import Crust

PLUME_RADIUS_KM = 800.
PLATEAU_KM = 30.        # oceanic plateau crust: about Ontong Java's
ACTIVE_LIP_MYR = 10.    # a LIP is 'active' up to this long after it erupts (its own step), 'former' after


def cap(lon: float, lat: float, radius_km: float) -> object:
    """A spherical cap of `radius_km` round (lon, lat)."""
    centre = sphere.from_lonlat(lon, lat)
    pole = sphere.unit(np.cross(centre, [0., 0., 1.]) if abs(centre[2]) < .99 else np.cross(centre, [1., 0., 0.]))
    start = sphere.rotate(centre, sphere.rotation(pole * radius_km / sphere.RADIUS_KM))
    ring = np.array([sphere.rotate(start, sphere.rotation(centre * 2. * math.pi * k / 72.)) for k in range(72)])
    return sphere.polygon(ring)


def plate_under(world, lon: float, lat: float):
    point = spherely.create_point(float(lon), float(lat))
    return next((p for p in world.plates() if spherely.contains(world.regions[p], point)), None)


def place_hotspot(hotspots: list, event: dict, time_myr: float) -> dict:
    ident = max((h['ident'] for h in hotspots), default=0) + 1
    hotspot = dict(ident=ident, name=event.get('name') or 'Hotspot %d' % ident, point=list(event['point']),
                   born_myr=time_myr, ended_myr=None, marks=[])
    hotspots.append(hotspot)
    return dict(kind='hotspot_placed', hotspot=ident, name=hotspot['name'])


def end_hotspot(hotspots: list, ident: int, time_myr: float) -> list:
    for h in hotspots:
        if h['ident'] == ident and h['ended_myr'] is None:
            h['ended_myr'] = time_myr
            return [dict(kind='hotspot_ended', hotspot=ident, name=h['name'])]
    return [dict(kind='declaration_missed', text='There is no active hotspot %d to end.' % ident)]


def emplace_lip(world, features: list, event: dict, time_myr: float) -> dict:
    """Erupt a LIP: plateau crust over ocean, a feature everywhere it lands; returns the event record."""
    from paradise.terranes import Feature, _next_ident
    lon, lat = event['point']
    province = cap(lon, lat, event.get('radius_km', PLUME_RADIUS_KM))
    name = event.get('name') or 'LIP'
    plateau, flood = 0., 0.
    for plate in world.plates():
        landed = sphere.intersection(province, world.regions[plate])
        if sphere.area_km2(landed) < sphere.DETAIL_KM2:
            continue
        local = sphere.rotate_geography(landed, world.frame(plate).T, tolerance=1e-4)
        made = []
        for c in world.crust:
            if c.plate != plate or c.kind != 'ocean':
                continue
            inside = sphere.intersection(c.shape, local)
            if sphere.area_km2(inside) >= sphere.DETAIL_KM2:
                c.shape = sphere.difference(c.shape, inside)
                made.append(inside)
        world.crust = [c for c in world.crust if sphere.area_km2(c.shape) >= sphere.DETAIL_KM2]
        if made:
            shape = sphere.union(*made)
            world.crust.append(Crust(plate, 'continent', shape, None, PLATEAU_KM))
            plateau += sphere.area_km2(shape)
        flood += sphere.area_km2(local) - sum(sphere.area_km2(m) for m in made)
        features.append(Feature(_next_ident(features), 'lip', plate, local, time_myr, name))
    return dict(kind='lip_erupted', name=name, plateau_km2=round(plateau), flood_km2=round(flood))


def mark_trails(world, hotspots: list, time_myr: float):
    """A trail mark for every active hotspot, on the plate over it now, in that plate's frame."""
    for h in hotspots:
        if h['ended_myr'] is not None:
            continue
        plate = plate_under(world, *h['point'])
        if plate is None:
            continue
        local = sphere.rotate(sphere.from_lonlat(*h['point']), world.frame(plate).T)
        h['marks'].append(dict(plate=plate, point=[float(v) for v in local], time_myr=time_myr))


def rehome_marks(world, hotspots: list):
    """Trail marks on a plate that has gone go to the plate now under them, as cratons do."""
    plates = set(world.plates())
    for h in hotspots:
        for m in h['marks']:
            if m['plate'] in plates:
                continue
            here = sphere.rotate(np.array(m['point']), world.frames.get(m['plate'], np.eye(3)))
            lon, lat = sphere.to_lonlat(here)
            new = plate_under(world, float(lon), float(lat))
            if new is not None:
                m['plate'], m['point'] = new, [float(v) for v in sphere.rotate(here, world.frame(new).T)]


def trails(hotspots: list) -> list:
    """[(hotspot, plate, [unit vectors in the plate's frame], first time, last time)]: each stretch of
    a trail on one plate, in order of time."""
    result = []
    for h in hotspots:
        run = []
        for m in h['marks'] + [None]:
            if run and (m is None or m['plate'] != run[-1]['plate']):
                if len(run) >= 2:
                    result.append((h, run[0]['plate'], np.array([r['point'] for r in run]),
                                   run[0]['time_myr'], run[-1]['time_myr']))
                run = []
            if m is not None:
                run.append(m)
    return result
