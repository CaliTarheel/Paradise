"""Elevation from the world's features, and the maps made from it.

This is a rendering, not a landscape evolution model: every term is a rule of thumb for what a
feature looks like, applied at the end of a step.

- **Ocean** deepens with age: 2.6 km at a ridge, then 0.345 km per square root of Myr, levelling
  out at 6.4 km (the plate-cooling depth-age curve).
- **Continent** stands from its thickness: 0.5 km above sea level at 35 km thick, 0.15 km higher
  per km thicker (Airy isostasy, as the core's ISOSTASY), and 0.4 km lower per km thinner, so the
  30 km crust of an oceanic plateau lies about 2 km down.
- **Uplift**:
  - active orogenies 3 km, former 1.5 km, old 0.5 km;
  - continental arcs 2.5 km;
  - island arcs stand at least ISLAND_ARC_KM, whatever their thin crust, so the noise makes an
    island chain of them.
- **Trenches**: 2.5 km deeper in a band TRENCH_KM wide on the downgoing side.

The base field is smoothed over SMOOTH_KM, which makes continental shelves and slopes out of the
steps between pieces. Uplift is smoothed separately and much less (UPLIFT_SMOOTH_KM): a 60 km arc
smoothed over 60 km would sink back into the sea. Fractal noise, rougher on high ground, then breaks up coasts and ridges.

The maps are a hypsometric colour image, a 16-bit heightmap, and Aerie's ROCKE-3D 72 by 46 grid,
both as its NPZ contract and as a TOPO NetCDF file in the form its studies use (all water is
ocean; no lakes or land ice).
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
import spherely
from PIL import Image
from scipy import ndimage

from paradise.core import sphere

RIDGE_KM, DEPTH_PER_ROOT_MYR, ABYSS_KM = 2.6, .345, 6.4
CONTINENT_BASE_KM, THICK_KM_PER_KM, THIN_KM_PER_KM, CONTINENT_KM = .5, .15, .4, 35.
# Every feature is read at a point displaced by fractal noise, up to about WARP_KM, so its edges
# (coasts above all) come out fractal rather than as the smooth outlines the model carries.
WARP_KM = 250.
UPLIFT_KM = dict(active=3., former=1.5, old=.5, continental_arc=2.5)
ISLAND_ARC_KM = .8           # island arcs stand at least this high, whatever their crust
UPLIFT_SMOOTH_KM = 15.       # mountain belts and arcs keep their width; only their edges soften
TRENCH_KM, TRENCH_DEEPEN_KM = 80., 2.5
SMOOTH_KM = 60.
# The core keeps the volume of shortened continent, and can pile it into crust far thicker than any
# real orogen. For the map, crust counts as no thicker than this, and land above SOFT_TOP_KM
# levels off towards SOFT_TOP_KM + SOFT_RANGE_KM (Tibet, not Everest, is the scale).
THICKEST_KM, SOFT_TOP_KM, SOFT_RANGE_KM = 75., 3., 3.
HEIGHT_RANGE_M = (-11000., 9000.)   # what the 16-bit heightmap spans


# ------------------------------------------------------------------------------------------- noise
def _lattice(ix, iy, iz, seed):
    h = (ix.astype(np.uint64) * np.uint64(73856093)) ^ (iy.astype(np.uint64) * np.uint64(19349663)) \
        ^ (iz.astype(np.uint64) * np.uint64(83492791)) ^ np.uint64(seed * 2654435761 & 0xffffffff)
    h ^= h >> np.uint64(13)
    h *= np.uint64(0x5bd1e995)
    h ^= h >> np.uint64(15)
    return (h & np.uint64(0xffff)).astype(np.float64) / 32767.5 - 1.


def value_noise(points, frequency, seed):
    """Smooth noise in [-1, 1] at unit vectors `points`, varying on a scale of 1/`frequency` radians."""
    p = np.asarray(points, float) * frequency + 1000.
    i = np.floor(p).astype(np.int64)
    f = p - i
    f = f * f * (3. - 2. * f)
    result = 0.
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = (f[..., 0] if dx else 1. - f[..., 0]) * (f[..., 1] if dy else 1. - f[..., 1]) \
                    * (f[..., 2] if dz else 1. - f[..., 2])
                result = result + w * _lattice(i[..., 0] + dx, i[..., 1] + dy, i[..., 2] + dz, seed)
    return result


def fbm(points, octaves=6, frequency=4., seed=1, ridged=False, persistence=.5):
    total, amplitude, norm = 0., 1., 0.
    for k in range(octaves):
        n = value_noise(points, frequency * 2. ** k, seed + 17 * k)
        if ridged:
            n = 1. - 2. * np.abs(n)
        total = total + amplitude * n
        norm += amplitude
        amplitude *= persistence
    return total / norm


# --------------------------------------------------------------------------------------- elevation
def _grid(width):
    height = width // 2
    lon = np.linspace(-180., 180., width, endpoint=False) + 180. / width
    lat = np.linspace(90., -90., height, endpoint=False) - 90. / height
    return np.meshgrid(lon, lat)


def _inside(geography, points, xyz=None):
    """Which of `points` (a spherely point array, made once and reused) lie in `geography`.

    With the points' unit vectors `xyz`, only those within a cap round the shape's vertices are
    tested: most shapes are small, and testing every point on a plate against each was most of the
    time a map took.
    """
    spherely.prepare(geography)
    if xyz is not None:
        rings = [r for poly in sphere.polygons(geography) for r in poly]
        if rings:
            vertices = np.vstack(rings)
            centre = vertices.sum(axis=0)
            if np.linalg.norm(centre) > 1e-6:
                centre = centre / np.linalg.norm(centre)
                reach = float(np.min(vertices @ centre))
                if reach > .2:   # the cap is well inside a hemisphere: worth it
                    near = xyz @ centre >= reach - 1e-6
                    result = np.zeros(len(xyz), bool)
                    if near.any():
                        result[near] = np.asarray(spherely.contains(geography, points[near]))
                    return result
    return np.asarray(spherely.contains(geography, points))


def _smooth(z, lat, km):
    """Gaussian smoothing over about `km` on the sphere: along meridians, then along each parallel."""
    height, width = z.shape
    km_per_px = math.pi * sphere.RADIUS_KM / height
    z = ndimage.gaussian_filter1d(z, km / km_per_px, axis=0, mode='nearest')
    out = np.empty_like(z)
    for row in range(height):
        sigma = km / (km_per_px * max(math.cos(math.radians(lat[row, 0])), .02))
        out[row] = ndimage.gaussian_filter1d(z[row], min(sigma, width / 4.), mode='wrap')
    return out


def elevation(snapshot, width=1440, seed=1) -> tuple:
    """(elevation in m on a `width` by `width`/2 equirectangular grid, lon, lat) for a snapshot."""
    from paradise.terranes import offset_line
    world = snapshot.world
    lon, lat = _grid(width)
    xyz = sphere.from_lonlat(lon.ravel(), lat.ravel())
    noise_at = xyz
    # Read the world at warped points: three noise fields displace each point by up to about WARP_KM.
    # Slow falloff (persistence 0.65) keeps the fine octaves strong enough to make a coast fractal.
    warp = np.stack([fbm(xyz, octaves=6, frequency=6., seed=seed + 7 * k, persistence=.65) for k in range(3)], axis=-1)
    xyz = sphere.unit(xyz + warp * (2. * WARP_KM / sphere.RADIUS_KM))
    flat_lon, flat_lat = sphere.to_lonlat(xyz)
    z = np.full(flat_lon.shape, -RIDGE_KM)
    lift = np.zeros(flat_lon.shape)
    everywhere = spherely.points(flat_lon, flat_lat)
    for plate, region in world.plate_polygons().items():
        inside = _inside(region, everywhere)
        if not inside.any():
            continue
        local = sphere.rotate(xyz[inside], world.frame(plate).T)
        llon, llat = sphere.to_lonlat(local)
        here_points = spherely.points(llon, llat)
        here = np.full(inside.sum(), -ABYSS_KM)
        for piece in (c for c in world.crust if c.plate == plate):
            hit = _inside(piece.shape, here_points, local)
            if not hit.any():
                continue
            if piece.kind == 'continent':
                mesh = getattr(piece, 'mesh', None)
                thick = np.asarray(mesh.sample(local[hit]) if mesh is not None else piece.thickness_km, float)
                thick = np.minimum(thick * np.ones(int(hit.sum())), THICKEST_KM)
                here[hit] = np.where(thick >= CONTINENT_KM,
                                     CONTINENT_BASE_KM + THICK_KM_PER_KM * (thick - CONTINENT_KM),
                                     CONTINENT_BASE_KM - THIN_KM_PER_KM * (CONTINENT_KM - thick))
            else:
                age = 0. if piece.born_myr is None else max(world.time_myr - piece.born_myr, 0.)
                here[hit] = -min(RIDGE_KM + DEPTH_PER_ROOT_MYR * math.sqrt(age), ABYSS_KM)
        uplift = np.zeros(inside.sum())
        floor = np.full(inside.sum(), -np.inf)
        for belt in getattr(snapshot, 'orogens', []):
            if belt['plate'] == plate:
                hit = _inside(belt['shape'], here_points, local)
                uplift[hit] = np.maximum(uplift[hit], UPLIFT_KM[belt['kind']])
        for f in getattr(snapshot, 'features', []):
            if f.plate == plate and f.kind in UPLIFT_KM:
                hit = _inside(f.shape, here_points, local)
                uplift[hit] = np.maximum(uplift[hit], UPLIFT_KM[f.kind])
            elif f.plate == plate and f.kind == 'island_arc':
                floor[_inside(f.shape, here_points, local)] = ISLAND_ARC_KM
        z[inside] = here
        lift[inside] = np.maximum(uplift, np.where(np.isfinite(floor), floor - here, 0.))
    # Trenches: a band on the downgoing side.
    for t in getattr(snapshot, 'trenches', []):
        points = sphere.from_lonlat(*np.array(t['points']).T)
        if len(points) < 2:
            continue
        on_left = t['downgoing'] == t['plates'][0]
        far = offset_line(points, TRENCH_KM, to_left=on_left)
        if len(far) != len(points):
            continue
        ring = np.vstack((points, far[::-1])) if not on_left else np.vstack((far, points[::-1]))
        try:
            band = sphere.polygon_robust(ring)
        except (ValueError, RuntimeError):
            continue
        if sphere.area_km2(band) > .5 * sphere.SPHERE_KM2:
            continue
        z[_inside(band, everywhere, xyz)] -= TRENCH_DEEPEN_KM
    z = _smooth(z.reshape(lon.shape), lat, SMOOTH_KM) + _smooth(lift.reshape(lon.shape), lat, UPLIFT_SMOOTH_KM)
    over = np.clip(z - SOFT_TOP_KM, 0., None)
    z = np.where(z > SOFT_TOP_KM, SOFT_TOP_KM + SOFT_RANGE_KM * (1. - np.exp(-over / SOFT_RANGE_KM)), z)
    # Noise: gentle on the abyssal plain, more at the coast, most on high ground.
    rough = fbm(noise_at, octaves=6, frequency=6., seed=seed).reshape(lon.shape)
    ridged = fbm(noise_at, octaves=5, frequency=10., seed=seed + 101, ridged=True).reshape(lon.shape)
    coast = np.exp(-(z / .3) ** 2)
    high = np.clip(z, 0., None)
    z = z + .15 * rough + .45 * coast * rough + .35 * high * (ridged - .3)
    return z * 1000., lon, lat


# -------------------------------------------------------------------------------------------- maps
TINTS = [(-11000., (8, 20, 60)), (-6000., (20, 45, 100)), (-3000., (40, 80, 140)), (-200., (90, 150, 200)),
         (0., (150, 200, 225)), (1., (90, 140, 80)), (300., (130, 165, 95)), (1000., (190, 175, 120)),
         (2000., (160, 120, 80)), (3500., (130, 110, 100)), (5000., (240, 240, 240))]


def colour(z) -> np.ndarray:
    levels = np.array([t for t, _ in TINTS])
    colours = np.array([c for _, c in TINTS], float)
    out = np.stack([np.interp(z, levels, colours[:, k]) for k in range(3)], axis=-1)
    # Hill shading from the north-west, stronger on land.
    gy, gx = np.gradient(z)
    shade = np.clip(1. - .0006 * (gx - gy), .55, 1.35)
    shade = np.where(z > 0., shade, 1. + .3 * (shade - 1.))
    return np.clip(out * shade[..., None], 0, 255).astype(np.uint8)


def aerie_grid(z, lon, lat) -> dict:
    """Aerie's 72 by 46 ROCKE-3D grid (its NPZ contract) from an equirectangular elevation field.

    Each target cell takes the area-weighted mean of the samples in it. The polar rows are one cap
    each (beyond 88 degrees), repeated at all 72 longitudes with a 72nd of its area, as ModelE does.
    """
    lon_bounds = np.linspace(-180., 180., 73)
    lat_bounds = np.r_[-90., np.arange(-88., 89., 4.), 90.]
    lon_centres = lon_bounds[:-1] + 2.5
    lat_centres = np.r_[-90., np.arange(-86., 87., 4.), 90.]
    weight = np.cos(np.radians(lat))
    j = np.clip(np.searchsorted(lat_bounds, lat, side='right') - 1, 0, 45)
    i = np.clip(((lon + 180.) // 5.).astype(int), 0, 71)
    polar = (j == 0) | (j == 45)
    i = np.where(polar, 0, i)
    land = z > 0.
    out = {k: np.full((46, 72), np.nan) for k in ('land_fraction', 'mean_signed_elevation_m', 'mean_land_elevation_m',
                                                  'mean_ocean_elevation_m', 'dominant_surface_elevation_m')}
    flat = j * 72 + i
    total = np.bincount(flat.ravel(), weight.ravel(), 46 * 72)
    land_w = np.bincount(flat.ravel(), (weight * land).ravel(), 46 * 72)
    signed = np.bincount(flat.ravel(), (weight * z).ravel(), 46 * 72)
    land_z = np.bincount(flat.ravel(), (weight * z * land).ravel(), 46 * 72)
    ocean_z = np.bincount(flat.ravel(), (weight * z * ~land).ravel(), 46 * 72)
    with np.errstate(invalid='ignore', divide='ignore'):
        lf = land_w / total
        grids = dict(land_fraction=lf, mean_signed_elevation_m=signed / total,
                     mean_land_elevation_m=np.where(land_w > 0, land_z / land_w, np.nan),
                     mean_ocean_elevation_m=np.where(total - land_w > 0, ocean_z / (total - land_w), np.nan))
    for key, values in grids.items():
        values = values.reshape(46, 72)
        values[0, :] = values[0, 0]     # one cap, repeated
        values[45, :] = values[45, 0]
        out[key] = values
    out['dominant_surface_elevation_m'] = np.where(out['land_fraction'] >= .5, out['mean_land_elevation_m'],
                                                   out['mean_ocean_elevation_m'])
    band = 2. * math.pi * (np.sin(np.radians(lat_bounds[1:])) - np.sin(np.radians(lat_bounds[:-1]))) / 72.
    out['cell_area_steradian'] = np.repeat(band[:, None], 72, axis=1)
    out.update(longitude_deg_e=lon_centres, latitude_deg_n=lat_centres, longitude_bounds_deg=lon_bounds,
               latitude_bounds_deg=lat_bounds)
    return out


def write_topo(grid: dict, path, title: str):
    """A ROCKE-3D TOPO NetCDF file in the form Aerie's studies use."""
    from scipy.io import netcdf_file
    lf = np.nan_to_num(grid['land_fraction'])
    land = np.nan_to_num(grid['mean_land_elevation_m'])
    ocean = np.nan_to_num(grid['mean_ocean_elevation_m'])
    focean = 1. - lf
    hocean = np.where(focean > 0., np.maximum(-ocean, 10.), 0.)
    zatmo = lf * np.maximum(land, 0.)
    with netcdf_file(path, 'w', version=1) as f:
        f.title = title
        f.source = 'Paradise'
        f.Conventions = 'CF-1.8'
        f.rocke3d_grid = 'Atm72x46 (4 degree latitude by 5 degree longitude)'
        f.remapping = 'area-weighted mean of an equirectangular elevation field'
        f.surface_assumption = 'all mapped water is ocean; no prescribed lakes or land ice'
        f.generator = 'paradise.relief'
        f.createDimension('lon', 72)
        f.createDimension('lat', 46)
        f.createDimension('bnds', 2)

        def variable(name, dims, values, kind='f', **attributes):
            v = f.createVariable(name, kind, dims)
            v[:] = values
            for k, a in attributes.items():
                setattr(v, k, a)

        lon_b, lat_b = grid['longitude_bounds_deg'], grid['latitude_bounds_deg']
        variable('lon', ('lon',), grid['longitude_deg_e'], 'd', standard_name='longitude',
                 long_name='longitude at grid-cell centre', units='degrees_east', axis='X', bounds='lon_bnds')
        variable('lat', ('lat',), grid['latitude_deg_n'], 'd', standard_name='latitude',
                 long_name='latitude at ModelE A-grid centre', units='degrees_north', axis='Y', bounds='lat_bnds')
        variable('lon_bnds', ('lon', 'bnds'), np.stack((lon_b[:-1], lon_b[1:]), axis=1), 'd', units='degrees_east')
        variable('lat_bnds', ('lat', 'bnds'), np.stack((lat_b[:-1], lat_b[1:]), axis=1), 'd', units='degrees_north')
        zero = np.zeros((46, 72))
        for name, values, long_name in (
                ('focean', focean, 'fraction of grid-cell area covered by ocean'),
                ('flake', zero, 'fraction of grid-cell area covered by lake'),
                ('fgrnd', lf, 'fraction of grid-cell area covered by non-glaciated land'),
                ('fgice', zero, 'fraction of grid-cell area covered by glaciated land')):
            variable(name, ('lat', 'lon'), values, long_name=long_name, units='1')
        for name, values, long_name in (
                ('zatmo', zatmo, 'mean atmospheric surface elevation above datum'),
                ('hocean', hocean, 'mean ocean depth over the ocean portion'),
                ('hlake', zero, 'mean lake depth over the lake portion'),
                ('hgice', zero, 'mean land-ice thickness over the glaciated portion'),
                ('zsolid', zatmo - focean * hocean, 'area-mean solid-ground elevation above datum')):
            variable(name, ('lat', 'lon'), values, long_name=long_name, units='m')


def write_maps(history, directory, width=1440, seed=1) -> list:
    """Every map of the latest step into `directory`; returns the files written."""
    os.makedirs(directory, exist_ok=True)
    snapshot = history.current
    age = history.project.age_ma(snapshot.time_myr)
    z, lon, lat = elevation(snapshot, width, seed)
    written = []
    path = os.path.join(directory, 'elevation.png')
    Image.fromarray(colour(z)).save(path)
    written.append(path)
    low, high = HEIGHT_RANGE_M
    path = os.path.join(directory, 'heightmap-16bit.png')
    Image.fromarray(np.round((np.clip(z, low, high) - low) / (high - low) * 65535.).astype(np.uint16)).save(path)
    written.append(path)
    grid = aerie_grid(z, lon, lat)
    path = os.path.join(directory, 'rocke3d-4x5.npz')
    np.savez(path, **grid)
    written.append(path)
    path = os.path.join(directory, 'rocke3d-TOPO.nc')
    write_topo(grid, path, '%s at %.0f Ma, from Paradise' % (history.project.name, age))
    written.append(path)
    land = float((np.cos(np.radians(lat)) * (z > 0)).sum() / np.cos(np.radians(lat)).sum())
    path = os.path.join(directory, 'maps.json')
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(dict(project=history.project.name, age_ma=age, width=width, seed=seed, land_share=land,
                       heightmap_range_m=list(HEIGHT_RANGE_M), highest_m=float(z.max()), deepest_m=float(z.min())),
                  f, indent=1)
    written.append(path)
    return written
