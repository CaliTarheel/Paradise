"""Which plate each node is on, and where the plates carry it over an interval."""
from __future__ import annotations

import numpy as np
import pygplates
from scipy.spatial import cKDTree


def plate_ids(project, age: float, mesh) -> np.ndarray:
    """The plate carrying each node at `age`: the plate ID of the crust under it.

    The crust features (continents, ocean strips, terranes) carry the plate IDs the rotation file
    moves. A hand-built project's plate outlines can be drawings only (Artifexia's are labelled
    2001, 2002... and never move), so they are not used. Land beats ocean; among ocean strips the
    youngest (the latest to form) wins. A node no crust covers takes its nearest covered node's.
    """
    from gplates_gospl.elevation import _inside, _polygons
    n = len(mesh.xyz)
    ids = np.full(n, -1, int)
    born = np.full(n, np.inf)
    for feature, polygon, centre, radius in _polygons(project, 'ocean', age):
        begin = feature.get_valid_time()[0]
        nodes, _ = _inside(mesh, polygon, centre, radius)
        younger = nodes[begin <= born[nodes]]
        born[younger] = begin
        ids[younger] = feature.get_reconstruction_plate_id()
    for role in ('continent', 'craton'):
        for feature, polygon, centre, radius in _polygons(project, role, age):
            nodes, _ = _inside(mesh, polygon, centre, radius)
            ids[nodes] = feature.get_reconstruction_plate_id()
    missing = ids < 0
    if missing.any() and (~missing).any():
        _, nearest = cKDTree(mesh.xyz[~missing]).query(mesh.xyz[missing])
        ids[missing] = ids[~missing][nearest]
    ids[ids < 0] = 0
    return ids


def carried(project, ids, nodes, older: float, younger: float) -> np.ndarray:
    """Where each node's plate carries it from `older` to `younger` (Ma): unit vectors."""
    out = np.array(nodes, float)
    for plate in np.unique(ids):
        rotation = project.rotation_model.get_rotation(younger, int(plate), older)
        matrix = _matrix(rotation)
        on = ids == plate
        out[on] = nodes[on] @ matrix.T
    return out / np.linalg.norm(out, axis=1, keepdims=True)


def _matrix(rotation) -> np.ndarray:
    if rotation.represents_identity_rotation():
        return np.eye(3)
    lat, lon, angle = rotation.get_lat_lon_euler_pole_and_angle_degrees()
    axis = np.asarray(pygplates.PointOnSphere(lat, lon).to_xyz(), float)
    a = np.radians(angle)
    k = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
    return np.eye(3) + np.sin(a) * k + (1. - np.cos(a)) * (k @ k)
