"""An elevation map of a reconstruction at one age, from its layers.

Relief is drawn from what each reconstructed feature is, by rules of thumb in the spirit of
Artifexian's hand-coloured maps, not by any physics:

- the ocean floor deepens with the age of its crust (Stein and Stein's GDH1), the age being how
  long ago the crust's strip formed (its begin time);
- continents stand a few hundred metres up, cratons higher, ringed by a shelf and a slope;
- orogenies rise towards their middles, highest while active, lower once former, lowest once old;
- island arcs rise out of the sea; LIPs build plateaus; trenches deepen the ocean beside them, and
  rifts lower the land along them.

goSPL then erodes and deposits; these maps only say where the land and the deep sea are, and how
high the ranges stand, at each stage.
"""
from __future__ import annotations

import warnings

import numpy as np
import pygplates
import spherely
from scipy.spatial import cKDTree

RADIUS_KM = 6371.
_KM_PER_S2_M = RADIUS_KM / spherely.EARTH_RADIUS_METERS

UNKNOWN_OCEAN_MYR = 80.          # crust no ocean layer covers is taken as this old
CONTINENT_M, CRATON_M, COAST_M = 300., 550., 40.
COAST_RAMP_KM = 250.             # land rises from the coast to its interior height over this
SHELF_M, SHELF_KM, SLOPE_KM = -120., 80., 220.
OROGENY = {'orogeny_active': (3500., 180.), 'orogeny_former': (1800., 200.), 'orogeny_old': (700., 220.)}
ARC_PEAK_M, ARC_FOOT_M, ARC_TAPER_KM = 900., -800., 60.
CONTINENTAL_ARC_M = 1200.
LIP = {'lip_active': (1200., -1500.), 'lip_former': (300., -2600.)}   # on land, plateau top at sea
TRENCH_M, TRENCH_KM = -7000., 90.
RIFT_M, RIFT_KM = -500., 80.


def ocean_depth_m(age_myr):
    """Sea floor depth (negative, m) for crust `age_myr` old: GDH1 (Stein and Stein 1992)."""
    age = np.clip(np.asarray(age_myr, float), 0., None)
    return -np.where(age < 20., 2600. + 365. * np.sqrt(age), 5651. - 2473. * np.exp(-.0278 * age))


def _taper(inside_km, width_km):
    """0 at a polygon's edge, rising smoothly to 1 at `width_km` inside it."""
    s = np.clip(np.asarray(inside_km) / width_km, 0., 1.)
    return s * s * (3. - 2. * s)


class Mesh:
    """The nodes as unit vectors, S2 points and a k-d tree, made once."""

    def __init__(self, nodes):
        self.xyz = np.asarray(nodes, float)
        lat = np.degrees(np.arcsin(np.clip(self.xyz[:, 2], -1., 1.)))
        lon = np.degrees(np.arctan2(self.xyz[:, 1], self.xyz[:, 0]))
        self.lon, self.lat = lon, lat
        self.points = spherely.points(lon, lat)
        self.tree = cKDTree(self.xyz)

    def near(self, centre, radius_km):
        """Indices of nodes within `radius_km` (great circle) of the unit vector `centre`."""
        chord = 2. * np.sin(min(radius_km / RADIUS_KM, np.pi) / 2.)
        return np.asarray(self.tree.query_ball_point(centre, chord + 1e-9), int)


def _ring(geometry):
    """(lon, lat) pairs of a pygplates polygon's exterior ring, closed duplicates removed."""
    latlon = np.asarray([q.to_lat_lon() for q in geometry.get_exterior_ring_points()], float)
    keep = np.r_[True, np.any(np.abs(np.diff(latlon, axis=0)) > 1e-9, axis=1)]
    latlon = latlon[keep]
    if len(latlon) > 1 and np.allclose(latlon[0], latlon[-1]):
        latlon = latlon[:-1]
    return [(float(lon), float(lat)) for lat, lon in latlon]


def _polygon(geometry):
    """An S2 polygon for a hand-drawn ring (either way round: the smaller side), or None."""
    ring = _ring(geometry)
    if len(ring) < 3:
        return None
    try:
        return spherely.create_polygon(ring, oriented=False)
    except (ValueError, RuntimeError):
        return None


def _polygons(project, role, age):
    """(feature, S2 polygon, centre unit vector, radius km) for each polygon of `role` at `age`."""
    out, skipped = [], 0
    for feature, geometry in project.reconstructed(role, age):
        if not isinstance(geometry, pygplates.PolygonOnSphere):
            continue
        polygon = _polygon(geometry)
        if polygon is None:
            skipped += 1
            continue
        centre = np.asarray(geometry.get_boundary_centroid().to_xyz(), float)
        ring = np.asarray([pygplates.PointOnSphere(lat, lon).to_xyz() for lon, lat in _ring(geometry)])
        radius = float(np.degrees(np.arccos(np.clip(ring @ centre, -1., 1.)).max()) * np.pi / 180. * RADIUS_KM)
        out.append((feature, polygon, centre, radius))
    if skipped:
        warnings.warn('%d %s polygons at %g Ma could not be read and were left out' % (skipped, role, age))
    return out


def _inside(mesh, polygon, centre, radius):
    """Indices of the nodes inside `polygon`, and their distance (km) in from its edge."""
    candidates = mesh.near(centre, radius + 1.)
    if not len(candidates):
        return candidates, np.zeros(0)
    hit = np.asarray(spherely.contains(polygon, mesh.points[candidates]), bool)
    nodes = candidates[hit]
    if not len(nodes):
        return nodes, np.zeros(0)
    edge = spherely.boundary(polygon)
    return nodes, np.asarray(spherely.distance(mesh.points[nodes], edge)) * _KM_PER_S2_M


def _lines(project, role, age):
    out = []
    for feature, geometry in project.reconstructed(role, age):
        if isinstance(geometry, pygplates.PolylineOnSphere):
            latlon = np.asarray(geometry.to_lat_lon_array(), float)
            if len(latlon) >= 2:
                out.append(spherely.create_linestring([(float(b), float(a)) for a, b in latlon]))
    return out


def _distance_to(mesh, geographies, within_km):
    """Each node's distance (km) to the nearest of `geographies`, or inf beyond `within_km`."""
    d = np.full(len(mesh.xyz), np.inf)
    for g in geographies:
        d = np.minimum(d, np.asarray(spherely.distance(mesh.points, g)) * _KM_PER_S2_M)
    d[d > within_km] = np.inf
    return d


def elevation(project, age: float, mesh: Mesh) -> dict:
    """{'z': elevation (m), 'land': continental crust (bool), 'crust_age': ocean crust age (Myr)}."""
    n = len(mesh.xyz)
    # Ocean floor by crust age: the youngest strip that covers a node formed it.
    born = np.full(n, np.inf)
    for feature, polygon, centre, radius in _polygons(project, 'ocean', age):
        begin = feature.get_valid_time()[0]
        if not np.isfinite(begin):
            continue
        nodes, _ = _inside(mesh, polygon, centre, radius)
        born[nodes] = np.minimum(born[nodes], begin)
    crust_age = np.where(np.isfinite(born), born - age, UNKNOWN_OCEAN_MYR)
    z = ocean_depth_m(crust_age)

    # Continents, then the interior rising from the coast.
    land = np.zeros(n, bool)
    for _, polygon, centre, radius in _polygons(project, 'continent', age):
        nodes, _ = _inside(mesh, polygon, centre, radius)
        land[nodes] = True
    craton = np.zeros(n, bool)
    for _, polygon, centre, radius in _polygons(project, 'craton', age):
        nodes, _ = _inside(mesh, polygon, centre, radius)
        craton[nodes] = True
    craton &= land
    if land.any() and (~land).any():
        sea_tree = cKDTree(mesh.xyz[~land])
        coast_km = sea_tree.query(mesh.xyz[land])[0] * RADIUS_KM
        interior = np.where(craton[land], CRATON_M, CONTINENT_M)
        z[land] = COAST_M + (interior - COAST_M) * _taper(coast_km, COAST_RAMP_KM)
        # The shelf and the slope down to the ocean floor.
        land_tree = cKDTree(mesh.xyz[land])
        off_km = land_tree.query(mesh.xyz[~land])[0] * RADIUS_KM
        sea = z[~land]
        slope = np.clip((off_km - SHELF_KM) / (SLOPE_KM - SHELF_KM), 0., 1.)
        z[~land] = np.where(off_km <= SHELF_KM, SHELF_M, SHELF_M + (sea - SHELF_M) * slope)
    elif land.all():
        z[:] = CONTINENT_M

    # Orogenies: added to the land, highest in the middle.
    for role, (height, width) in OROGENY.items():
        lift = np.zeros(n)
        for _, polygon, centre, radius in _polygons(project, role, age):
            nodes, inside_km = _inside(mesh, polygon, centre, radius)
            lift[nodes] = np.maximum(lift[nodes], height * _taper(inside_km + .5 * width / 3., width))
        z += np.where(land, lift, .5 * lift)

    # Island arcs out of the sea; on a continent, a volcanic range.
    for _, polygon, centre, radius in _polygons(project, 'arc', age):
        nodes, inside_km = _inside(mesh, polygon, centre, radius)
        s = _taper(inside_km + ARC_TAPER_KM / 3., ARC_TAPER_KM)
        sea = ~land[nodes]
        z[nodes[sea]] = np.maximum(z[nodes[sea]], ARC_FOOT_M + (ARC_PEAK_M - ARC_FOOT_M) * s[sea])
        z[nodes[~sea]] += CONTINENTAL_ARC_M * s[~sea]

    # LIPs: a lift on land, a plateau at sea.
    for role, (on_land, plateau) in LIP.items():
        for _, polygon, centre, radius in _polygons(project, role, age):
            nodes, _ = _inside(mesh, polygon, centre, radius)
            sea = ~land[nodes]
            z[nodes[~sea]] += on_land
            z[nodes[sea]] = np.maximum(z[nodes[sea]], plateau)

    # Trenches deepen the ocean beside them; rifts lower the land along them.
    trench = _distance_to(mesh, _lines(project, 'trench', age), TRENCH_KM)
    deep = np.isfinite(trench) & ~land
    z[deep] = np.minimum(z[deep], TRENCH_M + (z[deep] - TRENCH_M) * _taper(trench[deep], TRENCH_KM))
    rift = _distance_to(mesh, _lines(project, 'rift', age), RIFT_KM)
    low = np.isfinite(rift) & land
    z[low] += RIFT_M * (1. - _taper(rift[low], RIFT_KM))
    return dict(z=z, land=land, crust_age=crust_age)
