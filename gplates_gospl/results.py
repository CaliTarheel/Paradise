"""Read a goSPL run's output back onto the mesh it was built on.

Each MPI process writes its block (`h5/topology.pN.h5` coordinates, `h5/gospl.K.pN.h5` fields) in
its own node order, so blocks are matched to the mesh by coordinates. Where blocks overlap they
must agree; a node no block holds is an error.
"""
from __future__ import annotations

import glob
import os
import re

import h5py
import numpy as np
from scipy.spatial import cKDTree


def outputs(run_dir: str) -> list:
    """The output numbers a run has written."""
    found = {int(m.group(1)) for f in os.listdir(os.path.join(run_dir, 'h5'))
             if (m := re.match(r'gospl\.(\d+)\.p0\.h5$', f))}
    return sorted(found)


def read(run_dir: str, k: int, mesh_v, fields=('elev', 'erodep', 'FA')) -> dict:
    """{field: values on the mesh nodes} for output `k`; `mesh_v` is the mesh's `v` (m)."""
    tree = cKDTree(np.asarray(mesh_v, float))
    n = len(mesh_v)
    out = {f: np.full(n, np.nan) for f in fields}
    blocks = sorted(glob.glob(os.path.join(run_dir, 'h5', 'topology.p*.h5')))
    for block in blocks:
        rank = re.search(r'topology\.p(\d+)\.h5$', block).group(1)
        with h5py.File(block, 'r') as f:
            coords = f['coords'][:]
        distance, index = tree.query(coords)
        if distance.max() > 1e3:   # metres: a block node is not a mesh node
            raise ValueError('block %s of %s does not lie on the mesh' % (rank, run_dir))
        with h5py.File(os.path.join(run_dir, 'h5', 'gospl.%d.p%s.h5' % (k, rank)), 'r') as f:
            for field in fields:
                values = f[field][:, 0]
                seen = ~np.isnan(out[field][index])
                if seen.any() and not np.allclose(out[field][index][seen], values[seen], atol=1e-3, rtol=1e-6):
                    raise ValueError('blocks disagree on %s at output %d' % (field, k))
                out[field][index] = values
    for field, values in out.items():
        if np.isnan(values).any():
            raise ValueError('%d nodes have no %s at output %d' % (int(np.isnan(values).sum()), field, k))
    return out
