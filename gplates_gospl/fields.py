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
    """
    tree = tree if tree is not None else cKDTree(nodes)
    within = np.ones(len(nodes), bool) if within is None else np.asarray(within, bool)
    chord = 2. * np.sin(radius_km / 6371. / 2.)
    sigma2 = (chord / 2.) ** 2
    out = np.empty(len(nodes))
    for k, near in enumerate(tree.query_ball_point(nodes, chord)):
        near = np.asarray(near)
        near = near[within[near] == within[k]]
        w = np.exp(-((nodes[near] - nodes[k]) ** 2).sum(1) / (2. * sigma2))
        out[k] = (values[near] * w).sum() / w.sum()
    return out
