"""A closed icosphere: quasi-uniform triangles on the unit sphere, no seam and no polar fan."""
from __future__ import annotations

import numpy as np

# Mean edge length, km, by subdivision level: 6 is 120 km (40,962 nodes), 7 is 60 km (163,842).
EDGE_KM = {5: 240., 6: 120., 7: 60., 8: 30.}


def icosphere(level: int) -> tuple:
    """(nodes (N, 3) unit vectors, triangles (M, 3) zero-based) for `level` subdivisions."""
    t = (1. + 5. ** .5) / 2.
    nodes = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t],
                      [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    nodes /= np.linalg.norm(nodes, axis=1, keepdims=True)
    tris = np.array([[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4],
                     [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8],
                     [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]])
    for _ in range(level):
        edges = np.sort(np.concatenate((tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]])), axis=1)
        unique, inverse = np.unique(edges, axis=0, return_inverse=True)
        mids = nodes[unique[:, 0]] + nodes[unique[:, 1]]
        mids /= np.linalg.norm(mids, axis=1, keepdims=True)
        m = len(tris)
        a, b, c = (len(nodes) + inverse.reshape(-1)[k * m:(k + 1) * m] for k in range(3))
        nodes = np.vstack((nodes, mids))
        v0, v1, v2 = tris[:, 0], tris[:, 1], tris[:, 2]
        tris = np.concatenate((np.stack((v0, a, c), 1), np.stack((v1, b, a), 1),
                               np.stack((v2, c, b), 1), np.stack((a, b, c), 1)))
    return nodes, tris
