"""Pictures of maps and goSPL results: flat maps and globes, coloured by height.

    python -m gplates_gospl.render BUILD_DIR RUN_DIR OUT.png [--every K]

shows, for every K-th output of the run in RUN_DIR (an `output` folder of the build in BUILD_DIR),
the map it was forced towards, goSPL's landscape, and its erosion and deposition with the main
rivers.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree

HEIGHT_STOPS = [(-8000, (10, 20, 60)), (-5000, (25, 50, 110)), (-3000, (45, 90, 160)), (-200, (90, 150, 210)),
                (0, (150, 200, 230)), (1, (60, 120, 60)), (300, (110, 160, 80)), (800, (190, 180, 110)),
                (2000, (150, 110, 70)), (3500, (235, 235, 235))]


def colour(values, stops=HEIGHT_STOPS):
    zs = np.array([s[0] for s in stops], float)
    cs = np.array([s[1] for s in stops], float)
    return np.stack([np.interp(values, zs, cs[:, k]) for k in range(3)], -1).astype(np.uint8)


def erosion_colour(ed):
    """Erosion red, deposition blue, by metres (saturating at 1 km)."""
    s = np.clip(np.asarray(ed) / 1000., -1., 1.)
    white = np.full(s.shape + (3,), 245.)
    red, blue = np.array([180., 40., 30.]), np.array([30., 70., 180.])
    out = np.where(s[..., None] < 0, white + (red - white) * -s[..., None], white + (blue - white) * s[..., None])
    return out.astype(np.uint8)


class Raster:
    """Mesh values sampled onto a longitude-latitude grid, inverse-distance weighted."""

    def __init__(self, nodes, width=1440):
        nodes = np.asarray(nodes, float)
        nodes = nodes / np.linalg.norm(nodes, axis=1, keepdims=True)
        h = width // 2
        lon = np.radians(np.linspace(-180., 180., width, endpoint=False) + 180. / width)
        lat = np.radians(np.linspace(90., -90., h, endpoint=False) - 90. / h)
        lo, la = np.meshgrid(lon, lat)
        xyz = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1).reshape(-1, 3)
        d, self.index = cKDTree(nodes).query(xyz, k=3)
        self.weight = 1. / np.maximum(d, 1e-12)
        self.shape = (h, width)

    def __call__(self, values):
        values = np.asarray(values, float)
        return ((values[self.index] * self.weight).sum(1) / self.weight.sum(1)).reshape(self.shape)


def globe(image, lat0, lon0, size=420):
    """An orthographic view of an equirectangular image, centred on (lat0, lon0)."""
    H, W, _ = image.shape
    y, x = np.mgrid[-1:1:size * 1j, -1:1:size * 1j]
    r2 = x ** 2 + y ** 2
    inside = r2 <= 1
    z = np.sqrt(np.clip(1 - r2, 0, 1))
    la, lo = np.radians(lat0), np.radians(lon0)
    c = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    e = np.array([-np.sin(lo), np.cos(lo), 0.]) if abs(lat0) < 89 else np.array([0., 1., 0.])
    n = np.cross(c, e)
    p = z[..., None] * c + x[..., None] * e + (-y)[..., None] * n
    lat = np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1)))
    lon = np.degrees(np.arctan2(p[..., 1], p[..., 0]))
    px = ((lon + 180) / 360 * W).astype(int).clip(0, W - 1)
    py = ((90 - lat) / 180 * H).astype(int).clip(0, H - 1)
    out = np.full((size, size, 3), 255, np.uint8)
    out[inside] = image[py[inside], px[inside]]
    return Image.fromarray(out)


def rivers(image, flow, raster, top=.02):
    """Darken the cells carrying the most water (the top `top` share of flow accumulation on land)."""
    fa = raster(np.log10(np.maximum(flow, 1.)))
    wet = fa > 0.
    if not wet.any():
        return image
    cut = np.quantile(fa[wet], 1. - top)
    out = image.copy()
    out[wet & (fa >= cut)] = (20, 60, 160)
    return out


def sheet(build_dir, run_dir, path, every=1, width=1440):
    from gplates_gospl import results
    mesh = np.load(os.path.join(build_dir, 'input', 'mesh.npz'))
    maps = np.load(os.path.join(build_dir, 'input', 'maps.npz'))
    raster = Raster(mesh['v'], width)
    rows = []
    for k in results.outputs(run_dir)[::every]:
        got = results.read(run_dir, k, mesh['v'])
        land = got['elev'] > 0
        panels = [colour(raster(maps['z'][k])), rivers(colour(raster(got['elev'])), np.where(land, got['FA'], 0.), raster),
                  erosion_colour(raster(got['erodep']))]
        row = Image.new('RGB', (3 * width // 2, width // 4), (255, 255, 255))
        for j, p in enumerate(panels):
            row.paste(Image.fromarray(p).resize((width // 2, width // 4)), (j * width // 2, 0))
        ImageDraw.Draw(row).text((6, 4), '%g Ma: map, goSPL, erosion (red) and deposition (blue)' % maps['ages'][k],
                                 fill=(0, 0, 0))
        rows.append(row)
    out = Image.new('RGB', (rows[0].width, rows[0].height * len(rows)), (255, 255, 255))
    for k, r in enumerate(rows):
        out.paste(r, (0, k * r.height))
    out.save(path)
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('build')
    parser.add_argument('run')
    parser.add_argument('out')
    parser.add_argument('--every', type=int, default=1)
    args = parser.parse_args()
    print(sheet(args.build, args.run, args.out, args.every))


if __name__ == '__main__':
    main()
