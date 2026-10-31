"""Build a goSPL run from a GPlates project.

    python -m gplates_gospl.build PROJECT_DIR OUT_DIR --from 2000 --to 1150 --every 10 --level 6

writes, in OUT_DIR:

- `input/mesh.npz`: the icosphere (`v`, metres), its triangles (`c`) and the first stage's
  elevation (`z`), as goSPL's `npdata` reads them;
- `input/forcing_NNN.npz` for every interval: `hdisp`, the plates' displacement rate (m/yr, the
  chord each node travels over the interval), and `t`, the vertical rate that turns this stage's
  map, carried by the plates, into the next one's (m/yr);
- `input/rain.npz` (`r`, m/yr) and `input/maps.npz`: every stage's elevation, and the surface
  goSPL's `fast` mode should reach at each stage (`expected`), for comparison;
- `input.yml`, with erosion and deposition, and `forcing-check.yml`, goSPL's `fast` mode, which
  applies the forcing alone: its surface at each stage should be the map;

goSPL's time starts at 0 at the oldest stage and counts years forward: 2000 Ma is 0 and 1150 Ma
is 850,000,000.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from gplates_gospl import motion
from gplates_gospl.fields import advected, sample
from gplates_gospl.elevation import Mesh, elevation
from gplates_gospl.mesh import EDGE_KM, icosphere
from gplates_gospl.project import ARTIFEXIA, Project

RADIUS_M = 6371000.
RAIN_M_YR = 1.


def stage_ages(older: float, younger: float, every: float) -> list:
    """Stage boundaries from `older` to `younger` Ma, `every` Myr apart."""
    ages = [float(a) for a in np.arange(older, younger - 1e-9, -every)]
    if ages[-1] - younger > 1e-9:
        ages.append(float(younger))
    return ages


def drawn_at(age: float, younger: float) -> float:
    """The age a stage's map is drawn at: the youngest stage a hair older, so that the features
    that end exactly there (Artifexia's all end at 1150 Ma) are still in it."""
    return age + .01 if abs(age - younger) < 1e-9 else age


def build(project_dir, out, older=2000., younger=1150., every=10., level=6, dt=1e5, roles=ARTIFEXIA, log=print):
    project = Project(project_dir, roles)
    nodes, tris = icosphere(level)
    mesh = Mesh(nodes)
    ages = stage_ages(older, younger, every)
    os.makedirs(os.path.join(out, 'input'), exist_ok=True)
    maps, ids_at = [], []
    t0 = time.time()
    for age in ages:
        maps.append(elevation(project, drawn_at(age, younger), mesh)['z'])
        ids_at.append(motion.plate_ids(project, drawn_at(age, younger), mesh))
        log('map %7.2f Ma  land %4.1f%%  %.0f s' % (age, 100. * (maps[-1] > 0).mean(), time.time() - t0))
    np.savez_compressed(os.path.join(out, 'input', 'mesh'), v=nodes * RADIUS_M, c=tris, z=maps[0])
    np.savez_compressed(os.path.join(out, 'input', 'rain'), r=np.full(len(nodes), RAIN_M_YR))
    events, surface, expected = [], maps[0], [maps[0]]
    for k in range(len(ages) - 1):
        years = (ages[k] - ages[k + 1]) * 1e6
        moved = motion.carried(project, ids_at[k], nodes, ages[k], drawn_at(ages[k + 1], younger))
        hdisp = (moved - nodes) * RADIUS_M / years
        # goSPL lifts each node over the interval, then carries the surface with the plates at its
        # end: the lift is what the next map shows where the plate takes the node, less what is
        # there. "What is there" is goSPL's own surface, not the last map: carrying a surface
        # re-interpolates it, and a lift worked out from the map let that difference build up
        # stage on stage until, by 1600 Ma, the surface no longer resembled the maps at all.
        target = sample(mesh, maps[k + 1], moved)
        uplift = (target - surface) / years
        surface = advected(target, nodes, moved)
        expected.append(surface)
        name = 'forcing_%03d' % k
        np.savez_compressed(os.path.join(out, 'input', name), hdisp=hdisp, t=uplift)
        start, end = (older - ages[k]) * 1e6, (older - ages[k + 1]) * 1e6
        events.append(dict(start=start, end=end, hdisp=['input/' + name, 'hdisp'], upsub=['input/' + name, 't']))
    np.savez_compressed(os.path.join(out, 'input', 'maps'), ages=np.array(ages), z=np.array(maps),
                        expected=np.array(expected), plates=np.array(ids_at))
    total = (older - ages[-1]) * 1e6
    for fast, name, output in ((False, 'input.yml', 'output'), (True, 'forcing-check.yml', 'output-check')):
        config = dict(
            name='%s from %g to %g Ma' % (os.path.basename(os.path.normpath(project_dir)), older, ages[-1]),
            domain=dict(npdata=['input/mesh', 'v', 'c', 'z'], flowdir=3, fast=fast, advect='interp'),
            time=dict(start=0., end=total, tout=every * 1e6, dt=dt),
            spl=dict(K=1.e-6, d=0.42, m=0.5),
            diffusion=dict(hillslopeKa=0.5, hillslopeKm=2.0, nonlinKm=100., clinSlp=1.e-5),
            sea=dict(position=0.),
            tectonics=events,
            climate=[dict(start=0., map=['input/rain', 'r'])],
            output=dict(dir=output, makedir=False),
        )
        with open(os.path.join(out, name), 'w') as f:
            json.dump(config, f, indent=1)   # JSON is YAML
    with open(os.path.join(out, 'stages.json'), 'w') as f:
        json.dump(dict(project=os.path.abspath(project_dir), level=level, edge_km=EDGE_KM.get(level),
                       nodes=len(nodes), ages_ma=ages,
                       gospl_years=[(older - a) * 1e6 for a in ages]), f, indent=1)
    return ages, maps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('project')
    parser.add_argument('out')
    parser.add_argument('--from', dest='older', type=float, default=2000.)
    parser.add_argument('--to', dest='younger', type=float, default=1150.)
    parser.add_argument('--every', type=float, default=10.)
    parser.add_argument('--level', type=int, default=6)
    parser.add_argument('--dt', type=float, default=1e5)
    args = parser.parse_args()
    build(args.project, args.out, args.older, args.younger, args.every, args.level, args.dt)


if __name__ == '__main__':
    main()
