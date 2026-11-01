"""Rain shadows: a coarse climate's rainfall redistributed over the terrain beneath it.

ROCKE-3D's cells are about 450 km across, and goSPL's orographic precipitation works only on flat
(2D) meshes, so a global run would get no rain shadows at all. `orographic` reshapes the rainfall
over the current surface:

- wind blowing up a slope wrings out more rain (windward slopes wetter);
- land with higher ground upwind of it lies in that ground's shadow (drier);
- the climate's rainfall averaged over a ROCKE-3D cell's width is kept: rain is moved, never made.

The wind comes from `winds`: by default the planet's latitude bands (easterly trades within 30
degrees of the equator, westerlies to 60, polar easterlies beyond), until the climate supplies
its own surface winds.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

RADIUS_KM = 6371.
LIFT_PER_SLOPE = 40.        # rain multiplied by 1 + this x the windward slope (m/m): 1 km in 25 km triples it
SHADOW_M = 1500.            # rain falls by e for each this much higher ground upwind
REACH_KM, STEP_KM = 600., 50.


def winds(nodes):
    """Tangent vectors the way the wind blows, by latitude band: easterly trades within 30 degrees of
    the equator, westerlies to 60, polar easterlies beyond. The wind dies away towards each band's
    edge (the calms of the horse latitudes and the polar front), so the rain does not change at
    one line of latitude."""
    nodes = np.asarray(nodes, float)
    east = np.cross([0., 0., 1.], nodes)
    norm = np.linalg.norm(east, axis=1, keepdims=True)
    east = np.where(norm > 1e-9, east / np.maximum(norm, 1e-12), [1., 0., 0.])
    lat = np.abs(np.degrees(np.arcsin(np.clip(nodes[:, 2], -1., 1.))))
    # +1 westerly, -1 easterly, through 0 at 30 and 60 degrees
    blow = np.where(lat < 30., -np.clip((30. - lat) / 8., 0., 1.),
                    np.where(lat < 60., np.clip(np.minimum(lat - 30., 60. - lat) / 8., 0., 1.),
                             -np.clip((lat - 60.) / 8., 0., 1.)))
    return blow[:, None] * east


def _upwind(nodes, wind, km):
    """Each node's point `km` upwind, along the great circle the wind follows."""
    a = km / RADIUS_KM
    return np.cos(a) * nodes - np.sin(a) * wind


def orographic(rain, z, nodes, wind=None, keep_km=450.):
    """`rain` (m/yr, on the nodes) reshaped over the surface `z` (m).

    The rain is moved, not made: the result is rescaled so that its average over `keep_km` (a
    ROCKE-3D cell is about 450 km across) is the climate's. Rescaling each cell exactly left a seam
    at every cell edge.
    """
    from gplates_gospl.fields import smoothed
    nodes = np.asarray(nodes, float)
    nodes = nodes / np.linalg.norm(nodes, axis=1, keepdims=True)
    z = np.asarray(z, float)
    rain = np.asarray(rain, float)
    wind = winds(nodes) if wind is None else wind
    speed = np.linalg.norm(wind, axis=1)
    direction = wind / np.maximum(speed, 1e-12)[:, None]
    tree = cKDTree(nodes)

    def height(points):
        d, i = tree.query(points, k=3)
        w = 1. / np.maximum(d, 1e-12)
        return (z[i] * w).sum(1) / w.sum(1)

    ground = np.maximum(z, 0.)
    near = height(_upwind(nodes, direction, STEP_KM))
    slope = np.maximum(ground - np.maximum(near, 0.), 0.) / (STEP_KM * 1000.)
    barrier = np.zeros(len(z))
    for km in np.arange(STEP_KM, REACH_KM + 1e-9, STEP_KM):
        barrier = np.maximum(barrier, height(_upwind(nodes, direction, km)))
    shade = np.maximum(barrier - ground, 0.)
    factor = (1. + speed * LIFT_PER_SLOPE * slope) * np.exp(-speed * shade / SHADOW_M)
    shaped = rain * factor
    keep = smoothed(rain, nodes, keep_km) / np.maximum(smoothed(shaped, nodes, keep_km), 1e-12)
    return shaped * keep
