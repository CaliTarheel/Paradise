"""To and from ROCKE-3D's 72 x 46 grid, as Aerie's `map_definitions/rocke3d_4x5` defines it.

- `grid_npz` writes a world's elevation as Aerie's grid NPZ (land fraction, mean elevations,
  cell areas), the topography a ROCKE-3D world is set up from.
- `read_precipitation` reads the extractor's climate CSV (`j`, `i`, then `tNN_c`, `pNN_mm_day` for
  twelve months, then currents) into annual-mean rainfall on the grid, in m/yr, each month
  weighted by its length in days.
- `rain_on_mesh` interpolates that onto goSPL's mesh.
- `random_climate_csv` writes a made-up climate in the extractor's format, to prove the route
  before a real ROCKE-3D run exists.

Grid conventions (`grid-definition.json`): arrays are `[j, i]`, latitude south to north then
longitude west to east; `j = 0` and `j = 45` are the polar caps, one physical cell each,
repeated at all 72 longitudes. Longitude centres run -177.5 to 177.5 in 5 degree steps; latitude
centres are -90, -86 to 86 in 4 degree steps, and 90.
"""
from __future__ import annotations

import csv

import numpy as np

JM, IM = 46, 72
LON_CENTRES = -177.5 + 5. * np.arange(IM)
LON_BOUNDS = -180. + 5. * np.arange(IM + 1)
LAT_CENTRES = np.r_[-90., -86. + 4. * np.arange(JM - 2), 90.]
LAT_BOUNDS = np.r_[-90., -88. + 4. * np.arange(JM - 1), 90.]
SMAC_MONTH_DAYS = [44, 45, 44, 44, 45, 44, 44, 45, 44, 44, 45, 44]
DAYS_PER_YEAR = 365.25   # goSPL's year: rainfall in m/yr


def cell_of(nodes):
    """The grid cell (`j * 72 + i`) holding each unit vector; each polar cap is one cell."""
    nodes = np.asarray(nodes, float)
    nodes = nodes / np.linalg.norm(nodes, axis=1, keepdims=True)
    lat = np.degrees(np.arcsin(np.clip(nodes[:, 2], -1., 1.)))
    lon = np.degrees(np.arctan2(nodes[:, 1], nodes[:, 0]))
    j = np.clip(np.searchsorted(LAT_BOUNDS, lat, side='right') - 1, 0, JM - 1)
    i = np.clip(np.searchsorted(LON_BOUNDS, lon, side='right') - 1, 0, IM - 1)
    i = np.where((j == 0) | (j == JM - 1), 0, i)
    return j * IM + i


def cell_area_steradian():
    band = np.sin(np.radians(LAT_BOUNDS[1:])) - np.sin(np.radians(LAT_BOUNDS[:-1]))
    area = np.outer(band, np.full(IM, np.radians(5.)))
    area[0, :] = area[-1, :] = 2. * np.pi * band[0]   # each cap is one cell, repeated
    return area


def grid_npz(path, nodes, z, samples=8):
    """Write elevation `z` (m, on the unit vectors `nodes` of a mesh) as Aerie's grid NPZ.

    Each cell is sampled on a `samples` x `samples` grid of points, each taking its nearest node's
    elevation and weighted by its area; a cap is sampled round its whole circle. A mesh coarser
    than the grid's narrow polar cells (5 degrees of longitude is under 40 km at 86 degrees) still
    gives every cell a value.
    """
    from scipy.spatial import cKDTree
    nodes = np.asarray(nodes, float)
    nodes = nodes / np.linalg.norm(nodes, axis=1, keepdims=True)
    z = np.asarray(z, float)
    tree = cKDTree(nodes)
    f = (np.arange(samples) + .5) / samples
    shape = (JM, IM)
    out = {k: np.zeros(shape) for k in ('fraction', 'mean', 'land', 'sea')}
    for j in range(JM):
        lats = LAT_BOUNDS[j] + (LAT_BOUNDS[j + 1] - LAT_BOUNDS[j]) * f
        for i in range(IM if 0 < j < JM - 1 else 1):
            lons = (LON_BOUNDS[i] + 5. * f) if 0 < j < JM - 1 else (-180. + 360. * f)
            la, lo = np.meshgrid(np.radians(lats), np.radians(lons), indexing='ij')
            points = np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), -1).reshape(-1, 3)
            w = np.cos(la).reshape(-1)
            values = z[tree.query(points)[1]]
            land = values > 0.
            out['fraction'][j, i] = (w * land).sum() / w.sum()
            out['mean'][j, i] = (w * values).sum() / w.sum()
            out['land'][j, i] = (w * values)[land].sum() / w[land].sum() if land.any() else np.nan
            out['sea'][j, i] = (w * values)[~land].sum() / w[~land].sum() if (~land).any() else np.nan
    for a in out.values():
        a[0, :], a[-1, :] = a[0, 0], a[-1, 0]   # the caps, repeated
    dominant = np.where(out['fraction'] >= .5, out['land'], out['sea'])
    np.savez_compressed(path, longitude_deg_e=LON_CENTRES, latitude_deg_n=LAT_CENTRES,
                        longitude_bounds_deg=LON_BOUNDS, latitude_bounds_deg=LAT_BOUNDS,
                        cell_area_steradian=cell_area_steradian(), land_fraction=out['fraction'],
                        mean_signed_elevation_m=out['mean'], mean_land_elevation_m=out['land'],
                        mean_ocean_elevation_m=out['sea'], dominant_surface_elevation_m=dominant)


def read_precipitation(path, month_days=SMAC_MONTH_DAYS):
    """Annual-mean rainfall (m/yr) on the grid, `[j, i]`, from the extractor's climate CSV."""
    rain = np.full((12, JM, IM), np.nan)
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    if len(rows) != JM * IM:
        raise ValueError('%s has %d rows, not %d' % (path, len(rows), JM * IM))
    for row in rows:
        j, i = int(row['j']), int(row['i'])
        for m in range(12):
            rain[m, j, i] = float(row['p%02d_mm_day' % (m + 1)])
    if np.isnan(rain).any() or (rain < 0).any():
        raise ValueError('%s has missing or negative precipitation' % path)
    weights = np.asarray(month_days, float)[:, None, None]
    mm_day = (rain * weights).sum(0) / weights.sum()
    return mm_day * DAYS_PER_YEAR / 1000.


def rain_on_mesh(grid, nodes):
    """Grid values (`[j, i]`) at the unit vectors `nodes`: bilinear between cell centres, periodic in
    longitude, the caps taken at the poles."""
    nodes = np.asarray(nodes, float)
    nodes = nodes / np.linalg.norm(nodes, axis=1, keepdims=True)
    lat = np.degrees(np.arcsin(np.clip(nodes[:, 2], -1., 1.)))
    lon = np.degrees(np.arctan2(nodes[:, 1], nodes[:, 0]))
    j0 = np.clip(np.searchsorted(LAT_CENTRES, lat, side='right') - 1, 0, JM - 2)
    fj = np.clip((lat - LAT_CENTRES[j0]) / (LAT_CENTRES[j0 + 1] - LAT_CENTRES[j0]), 0., 1.)
    x = (lon - LON_CENTRES[0]) / 5.
    i0 = np.floor(x).astype(int) % IM
    fi = x - np.floor(x)
    i1 = (i0 + 1) % IM
    g = np.asarray(grid, float)
    low = g[j0, i0] * (1. - fi) + g[j0, i1] * fi
    high = g[j0 + 1, i0] * (1. - fi) + g[j0 + 1, i1] * fi
    return low * (1. - fj) + high * fj


def random_climate_csv(path, land_fraction=None, seed=0):
    """A made-up climate in the extractor's CSV format: a wet equator, dry subtropics and wetter
    mid-latitudes, broken up by random weather cells, wetter over sea, varying through the year.
    Temperatures are a plain latitude profile; currents are NaN."""
    rng = np.random.default_rng(seed)
    lat, lon = np.meshgrid(np.radians(LAT_CENTRES), np.radians(LON_CENTRES), indexing='ij')
    zonal = 1.2 + 3.5 * np.exp(-(np.degrees(lat) / 12.) ** 2) - 1.0 * np.exp(-((np.abs(np.degrees(lat)) - 25.) / 8.) ** 2) \
        + 0.8 * np.exp(-((np.abs(np.degrees(lat)) - 50.) / 12.) ** 2)
    xyz = np.stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)), -1)
    cells = np.zeros(lat.shape)
    for _ in range(40):
        centre = rng.normal(size=3)
        centre /= np.linalg.norm(centre)
        width = rng.uniform(.15, .45)
        cells += rng.normal() * np.exp(-((xyz - centre) ** 2).sum(-1) / (2 * width ** 2))
    annual = np.clip(zonal * np.exp(.5 * cells / max(np.abs(cells).max(), 1e-9)), .05, None)
    if land_fraction is not None:
        annual *= 1. - .35 * np.asarray(land_fraction)
    temperature = 30. * np.cos(lat) - 10.
    with open(path, 'w', newline='') as f:
        names = ['j', 'i'] + [n for m in range(1, 13) for n in ('t%02d_c' % m, 'p%02d_mm_day' % m)] + \
                ['u_face_cm_s', 'v_face_cm_s', 'w_top_cm_s']
        writer = csv.writer(f)
        writer.writerow(names)
        for j in range(JM):
            for i in range(IM):
                jj, ii = j, (0 if j in (0, JM - 1) else i)
                row = [j, i]
                for m in range(12):
                    season = 1. + .4 * np.sin(2 * np.pi * m / 12.) * np.sign(LAT_CENTRES[j])
                    row += ['%.4f' % temperature[jj, ii], '%.6f' % (annual[jj, ii] * season)]
                writer.writerow(row + ['nan', 'nan', 'nan'])
