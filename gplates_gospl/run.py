"""Run goSPL on a build interval by interval, keeping the landscape on the reconstruction.

    /opt/gospl-env/bin/python -m gplates_gospl.run BUILD_DIR RUN_DIR [--first K] [--count N] [--np 4]

Run in one go, goSPL erodes the land and the forcing never gives it back: on Artifexia the land
fell from 26% to 21% of the sphere, and its highest point from 1,490 m to 150 m, in the first
10 Myr. So each interval is run on its own, from where the last one ended, and its vertical forcing
is the build's plus the erosion and deposition so far, smoothed: what goSPL's surface lacks
against the surface the forcing alone would have reached (`expected`), averaged over
`--smooth-km`, land and sea apart. The regional height goes back to the reconstruction's; the
valleys, basins and drainage goSPL has cut stay. The erosion and deposition still to come over the
interval, going by the last interval's rates, are made up for as they happen; the first interval is
run once beforehand to measure them (`spin_up/`).

RUN_DIR keeps, per interval, goSPL's input and output (`interval_NNN/`), and the surface it ended
on (`state_NNN.npz`: elevation `z`, cumulative erosion and deposition `ed`, and the standing
`change` against the forcing alone, on the mesh).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

from gplates_gospl import results
from gplates_gospl import rain as rain_shadows
from gplates_gospl.fields import advected, sample, smoothed

RADIUS_M = 6371000.
SEDIMENT_SINK = 1. - 2300. / 3300.   # of a sediment load's thickness, the plate sinks this much


class _Tree:
    """The mesh as `fields.sample` takes it: anything with a k-d tree of its unit nodes."""

    def __init__(self, tree):
        self.tree = tree


def _interval(build_dir, run_dir, template, k, state, uplift, hdisp, cells, v, processes, dt, name, rain=None):
    """Run goSPL over interval `k` from `state` under `uplift`; its last output on the mesh."""
    event = template['tectonics'][k]
    years = event['end'] - event['start']
    folder = os.path.join(run_dir, name)
    os.makedirs(os.path.join(folder, 'input'), exist_ok=True)
    np.savez_compressed(os.path.join(folder, 'input', 'mesh'), v=v, c=cells, z=state)
    np.savez_compressed(os.path.join(folder, 'input', 'forcing'), hdisp=hdisp, t=uplift)
    np.savez_compressed(os.path.join(folder, 'input', 'rain'),
                        r=rain if rain is not None else np.load(os.path.join(build_dir, 'input', 'rain.npz'))['r'])
    config = dict(template)
    config['time'] = dict(template['time'], start=event['start'], end=event['end'], tout=years)
    if dt:
        config['time']['dt'] = dt
    config['tectonics'] = [dict(start=event['start'], end=event['end'], hdisp=['input/forcing', 'hdisp'],
                                upsub=['input/forcing', 't'])]
    config['climate'] = [dict(start=event['start'], map=['input/rain', 'r'])]
    config['output'] = dict(dir='output', makedir=False)
    with open(os.path.join(folder, 'input.yml'), 'w') as f:
        json.dump(config, f, indent=1)
    # Ubuntu's PETSc has no ParMETIS to split the mesh between processes; PETSc's own simple
    # partitioner stands in, unless the caller's PETSC_OPTIONS already say how.
    env = dict(os.environ)
    env.setdefault('PETSC_OPTIONS', '-petscpartitioner_type simple')
    command = ['mpirun', '--allow-run-as-root', '-x', 'PETSC_OPTIONS', '-np', str(processes),
               os.path.join(os.path.dirname(sys.executable), 'gospl'), '-i', 'input.yml']
    with open(os.path.join(folder, 'gospl.log'), 'w') as out:
        done = subprocess.run(command, cwd=folder, env=env, stdout=out, stderr=subprocess.STDOUT)
    if done.returncode:
        raise RuntimeError('goSPL failed on interval %d; see %s' % (k, os.path.join(folder, 'gospl.log')))
    output = os.path.join(folder, 'output')
    return results.read(output, results.outputs(output)[-1], v)


def run(build_dir, run_dir, first=0, count=None, processes=4, smooth_km=400., dt=None, erodibility=None, rain=None,
        orographic=False, log=print):
    mesh = np.load(os.path.join(build_dir, 'input', 'mesh.npz'))
    maps = np.load(os.path.join(build_dir, 'input', 'maps.npz'))
    with open(os.path.join(build_dir, 'input.yml')) as f:
        template = json.load(f)
    if erodibility is not None:
        template['spl'] = dict(template['spl'], K=erodibility)
    v, cells = mesh['v'], mesh['c']
    nodes = v / np.linalg.norm(v, axis=1, keepdims=True)
    tree = cKDTree(nodes)
    ages, expected = maps['ages'], maps['expected']
    base_rain = rain if rain is not None else np.load(os.path.join(build_dir, 'input', 'rain.npz'))['r']
    last = len(ages) - 1 if count is None else min(len(ages) - 1, first + count)
    os.makedirs(run_dir, exist_ok=True)
    if first == 0:
        state, ed, recent = mesh['z'], np.zeros(len(v)), None
    else:
        saved = np.load(os.path.join(run_dir, 'state_%03d.npz' % first))
        state, ed, recent = saved['z'], saved['ed'], saved['recent']
    for k in range(first, last):
        t0 = time.time()
        event = template['tectonics'][k]
        years = event['end'] - event['start']
        forcing = np.load(os.path.join(build_dir, event['upsub'][0] + '.npz'))
        hdisp = forcing['hdisp']
        land = expected[k] > 0.
        lacking = smoothed(expected[k] - state, nodes, smooth_km, within=land, tree=tree)
        if recent is None:
            here_rain = rain_shadows.orographic(base_rain, state, nodes) if orographic else base_rain
            # Spin-up: the first interval once as it is, for an erosion rate to anticipate.
            got = _interval(build_dir, run_dir, template, k, state, forcing['t'] + lacking / years, hdisp,
                            cells, v, processes, dt, 'spin_up', here_rain)
            # Its rates are where the plates took them; read back where the interval started.
            moved = v + hdisp * years
            recent = sample(_Tree(tree), got['erodep'] / years, moved / np.linalg.norm(moved, axis=1, keepdims=True))
        # What erosion and deposition will do over the interval, going by the last one, is made
        # up for as it happens: with only what they had done made up, every interval ended one
        # interval's erosion down, and the land shrank from 26% of the sphere to 21% in each.
        # On land, the erosion to come is made up for as it happens (smoothed, so the valleys stay).
        # At sea, sediment is not given room ahead of time: its weight presses the plate down by
        # SEDIMENT_SINK of its thickness (Airy isostasy, sediment 2,300 over mantle 3,300 kg/m3),
        # node by node, over the next interval. Smoothed like the land, deposits were let pile into
        # deep water: 0.7% of the sphere shoaled from below -1,000 m to above -200 m in 50 Myr.
        coming = np.where(land, smoothed(recent, nodes, smooth_km, within=land, tree=tree), 0.)
        sinking = np.where(land, 0., SEDIMENT_SINK * np.maximum(recent, 0.))
        uplift = forcing['t'] + lacking / years - coming - sinking
        # The corrections are smoothed, so they give a range back its average loss: crests, which
        # lose least, gained more than they lost, and Artifexia's highest point climbed 450 m an
        # interval past the map's. No node is lifted past the map's height where its plate takes
        # it (`expected` carried by the build's forcing); erosion may still take it below.
        target = expected[k] + forcing['t'] * years
        here_rain = base_rain
        if orographic:
            here_rain = rain_shadows.orographic(base_rain, state, nodes)
        uplift = np.minimum(uplift, np.maximum((target - state) / years, 0.))
        got = _interval(build_dir, run_dir, template, k, state, uplift, hdisp, cells, v, processes, dt,
                        'interval_%03d' % k, here_rain)
        moved = v + hdisp * years
        ed = advected(ed, v, moved) + got['erodep']
        recent = got['erodep'] / years
        state = got['elev']
        # `ed` adds up every interval's erosion and deposition, so where the forcing re-deepens a
        # trench each interval and sediment refills it, it counts every refill (15 km of it at
        # Artifexia's trenches). `change` is what stands: goSPL's surface against the surface the
        # forcing alone reaches, cut valleys negative, sediment cover positive.
        np.savez_compressed(os.path.join(run_dir, 'state_%03d' % (k + 1)), z=state, ed=ed, recent=recent,
                            change=state - expected[k + 1])
        d = state - maps['z'][k + 1]
        log('%7.2f Ma  %4.0f s  land %4.1f%% (map %4.1f%%)  top %5.0f m (map %5.0f)  rms from map %4.0f m'
            '  eroded on land: median %3.0f m'
            % (ages[k + 1], time.time() - t0, 100 * (state > 0).mean(), 100 * (maps['z'][k + 1] > 0).mean(),
               state.max(), maps['z'][k + 1].max(), np.sqrt((d ** 2).mean()),
               np.median(-got['erodep'][state > 0]) if (state > 0).any() else 0.))
    return state


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('build')
    parser.add_argument('run')
    parser.add_argument('--first', type=int, default=0)
    parser.add_argument('--count', type=int)
    parser.add_argument('--np', dest='processes', type=int, default=4)
    parser.add_argument('--smooth-km', type=float, default=400.)
    parser.add_argument('--dt', type=float)
    parser.add_argument('--rain', help="rainfall on the mesh (m/yr): an .npz with 'r', or a ROCKE-3D climate CSV "
                                       "(Aerie's extractor format), interpolated onto the mesh")
    parser.add_argument('--orographic', action='store_true',
                        help='reshape the rain over the terrain each interval: rain shadows (gplates_gospl.rain)')
    parser.add_argument('--K', dest='erodibility', type=float,
                        help="goSPL's stream-power erodibility; the build's is 1e-6, set for a 30 km mesh")
    args = parser.parse_args()
    rain = None
    if args.rain:
        v = np.load(os.path.join(args.build, 'input', 'mesh.npz'))['v']
        if args.rain.endswith('.npz'):
            rain = np.load(args.rain)['r']
        else:
            from gplates_gospl import rocke3d
            rain = rocke3d.rain_on_mesh(rocke3d.read_precipitation(args.rain), v)
    run(args.build, args.run, args.first, args.count, args.processes, args.smooth_km, args.dt, args.erodibility, rain,
        args.orographic,
        log=lambda s: print(s, flush=True))


if __name__ == '__main__':
    main()
