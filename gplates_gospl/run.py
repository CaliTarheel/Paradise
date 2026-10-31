"""Run goSPL on a build interval by interval, keeping the landscape on the reconstruction.

    /opt/gospl-env/bin/python -m gplates_gospl.run BUILD_DIR RUN_DIR [--first K] [--count N] [--np 4]

Run in one go, goSPL erodes the land and the forcing never gives it back: on Artifexia the land
fell from 26% to 21% of the sphere, and its highest point from 1,490 m to 150 m, in the first
10 Myr. So each interval is run on its own, from where the last one ended, and its vertical forcing
is the build's plus the erosion and deposition so far, smoothed: what goSPL's surface lacks
against the surface the forcing alone would have reached (`expected`), averaged over
`--smooth-km`, land and sea apart. The regional height goes back to the reconstruction's; the
valleys, basins and drainage goSPL has cut stay.

RUN_DIR keeps, per interval, goSPL's input and output (`interval_NNN/`), and the surface it ended
on (`state_NNN.npz`: elevation `z` and cumulative erosion and deposition `ed` on the mesh).
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
from gplates_gospl.fields import advected, smoothed

RADIUS_M = 6371000.


def run(build_dir, run_dir, first=0, count=None, processes=4, smooth_km=400., dt=None, log=print):
    mesh = np.load(os.path.join(build_dir, 'input', 'mesh.npz'))
    maps = np.load(os.path.join(build_dir, 'input', 'maps.npz'))
    with open(os.path.join(build_dir, 'input.yml')) as f:
        template = json.load(f)
    v, cells = mesh['v'], mesh['c']
    nodes = v / np.linalg.norm(v, axis=1, keepdims=True)
    tree = cKDTree(nodes)
    ages, expected = maps['ages'], maps['expected']
    last = len(ages) - 1 if count is None else min(len(ages) - 1, first + count)
    os.makedirs(run_dir, exist_ok=True)
    if first == 0:
        state, ed = mesh['z'], np.zeros(len(v))
    else:
        saved = np.load(os.path.join(run_dir, 'state_%03d.npz' % first))
        state, ed = saved['z'], saved['ed']
    for k in range(first, last):
        t0 = time.time()
        event = template['tectonics'][k]
        years = event['end'] - event['start']
        forcing = np.load(os.path.join(build_dir, event['upsub'][0] + '.npz'))
        hdisp = forcing['hdisp']
        land = expected[k] > 0.
        lacking = smoothed(expected[k] - state, nodes, smooth_km, within=land, tree=tree)
        uplift = forcing['t'] + lacking / years
        folder = os.path.join(run_dir, 'interval_%03d' % k)
        os.makedirs(os.path.join(folder, 'input'), exist_ok=True)
        np.savez_compressed(os.path.join(folder, 'input', 'mesh'), v=v, c=cells, z=state)
        np.savez_compressed(os.path.join(folder, 'input', 'forcing'), hdisp=hdisp, t=uplift)
        np.savez_compressed(os.path.join(folder, 'input', 'rain'),
                            r=np.load(os.path.join(build_dir, 'input', 'rain.npz'))['r'])
        config = dict(template)
        config['name'] = '%s, %g to %g Ma' % (template['name'], ages[k], ages[k + 1])
        config['time'] = dict(template['time'], start=event['start'], end=event['end'], tout=years)
        if dt:
            config['time']['dt'] = dt
        config['tectonics'] = [dict(start=event['start'], end=event['end'], hdisp=['input/forcing', 'hdisp'],
                                    upsub=['input/forcing', 't'])]
        config['climate'] = [dict(start=event['start'], map=['input/rain', 'r'])]
        config['output'] = dict(dir='output', makedir=False)
        with open(os.path.join(folder, 'input.yml'), 'w') as f:
            json.dump(config, f, indent=1)
        env = dict(os.environ, PETSC_OPTIONS='-petscpartitioner_type simple')
        command = ['mpirun', '--allow-run-as-root', '-x', 'PETSC_OPTIONS', '-np', str(processes),
                   os.path.join(os.path.dirname(sys.executable), 'gospl'),
                   '-i', 'input.yml']
        with open(os.path.join(folder, 'gospl.log'), 'w') as out:
            done = subprocess.run(command, cwd=folder, env=env, stdout=out, stderr=subprocess.STDOUT)
        if done.returncode:
            raise RuntimeError('goSPL failed on interval %d; see %s' % (k, os.path.join(folder, 'gospl.log')))
        got = results.read(os.path.join(folder, 'output'), results.outputs(os.path.join(folder, 'output'))[-1], v)
        moved = v + hdisp * years
        ed = advected(ed, v, moved) + got['erodep']
        state = got['elev']
        np.savez_compressed(os.path.join(run_dir, 'state_%03d' % (k + 1)), z=state, ed=ed)
        d = state - maps['z'][k + 1]
        log('%7.2f Ma  %4.0f s  land %4.1f%% (map %4.1f%%)  top %5.0f m (map %5.0f)  rms from map %4.0f m'
            % (ages[k + 1], time.time() - t0, 100 * (state > 0).mean(), 100 * (maps['z'][k + 1] > 0).mean(),
               state.max(), maps['z'][k + 1].max(), np.sqrt((d ** 2).mean())))
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
    args = parser.parse_args()
    run(args.build, args.run, args.first, args.count, args.processes, args.smooth_km, args.dt,
        log=lambda s: print(s, flush=True))


if __name__ == '__main__':
    main()
