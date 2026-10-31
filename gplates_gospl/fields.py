"""Fields on the mesh: sampling them elsewhere, carrying them as goSPL does, and smoothing them.

Nothing here needs pygplates, so goSPL's own environment can import it.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def sample(mesh, values, points, k=3):
    """`values` (on the mesh nodes) at unit vectors `points`, inverse-distance weighted; `mesh` is
    anything with a k-d tree of the unit node vectors as `tree`."""
    d, i = mesh.tree.query(points, k=k)
    w = 1. / np.maximum(d, 1e-12)
    return (values[i] * w).sum(1) / w.sum(1)


def advected(values, nodes, moved):
    """`values` at the `moved` nodes, interpolated back onto `nodes` as goSPL's `interp` scheme does
    it (`Tectonics._readAdvectionData`, `_advectPlates`): the three nearest moved nodes, weighted by
    inverse distance, or the one that lands exactly on a node."""
    d, i = cKDTree(moved).query(nodes, k=3)
    w = np.divide(1., d, out=np.zeros_like(d), where=d != 0)
    out = (w * values[i]).sum(1) / np.where(w.sum(1) > 0, w.sum(1), 1.)
    exact = d[:, 0] == 0
    out[exact] = values[i[exact, 0]]
    return out


def smoothed(values, nodes, radius_km, within=None, tree=None):
    """`values` averaged over a cap of `radius_km` round each node, weighted by a Gaussian.

    With `within` (a boolean per node), each node averages only the nodes of its own kind, so land
    and sea, say, are smoothed apart and a coast is not blurred across.

    The averaging is done on a coarse helper mesh, a few nodes to the radius: each node's value is
    pooled into its nearest helper node, the pools are smoothed (normalised convolution, so a pool
    of one kind is not diluted by the other's empty helper nodes), and the result is read back at
    each node. A fine mesh's own neighbourhoods (thousands of nodes each at 15 km) made the direct
    sum unworkable; this costs the same at any resolution. `tree` is accepted and ignored.
    """
    from gplates_gospl.mesh import icosphere
    nodes = np.asarray(nodes, float)
    nodes = nodes / np.linalg.norm(nodes, axis=1, keepdims=True)
    values = np.asarray(values, float)
    within = np.ones(len(nodes), bool) if within is None else np.asarray(within, bool)
    level = _helper_level(radius_km)
    helper, _ = icosphere(level)
    pool = cKDTree(helper).query(nodes)[1]
    weights = _gaussian(helper, radius_km)
    out = np.empty(len(nodes))
    for kind in (True, False):
        mine = within == kind
        if not mine.any():
            continue
        total = np.bincount(pool[mine], values[mine], minlength=len(helper))
        count = np.bincount(pool[mine], minlength=len(helper)).astype(float)
        total, count = weights @ total, weights @ count
        field = np.divide(total, count, out=np.full(len(helper), np.nan), where=count > 1e-9)
        known = np.isfinite(field)
        d, i = cKDTree(helper[known]).query(nodes[mine], k=3)
        w = 1. / np.maximum(d, 1e-12)
        out[mine] = (field[known][i] * w).sum(1) / w.sum(1)
    return out


def _helper_level(radius_km):
    """The coarsest icosphere with at least three nodes across the radius."""
    from gplates_gospl.mesh import EDGE_KM
    for level in sorted(EDGE_KM):
        if EDGE_KM[level] <= radius_km / 3.:
            return level
    return max(EDGE_KM)


def _gaussian(nodes, radius_km):
    """The sparse matrix of Gaussian weights (sigma half the radius) between nodes within `radius_km`."""
    from scipy import sparse
    chord = 2. * np.sin(radius_km / 6371. / 2.)
    pairs = cKDTree(nodes).query_pairs(chord, output_type='ndarray')
    rows = np.r_[pairs[:, 0], pairs[:, 1], np.arange(len(nodes))]
    cols = np.r_[pairs[:, 1], pairs[:, 0], np.arange(len(nodes))]
    d2 = ((nodes[rows] - nodes[cols]) ** 2).sum(1)
    return sparse.csr_matrix((np.exp(-d2 / (2. * (chord / 2.) ** 2)), (rows, cols)), shape=(len(nodes),) * 2)
