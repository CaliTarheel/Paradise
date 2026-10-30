"""Rigid continents meeting: weld two continental plates before their continents overlap.

A hand-built reconstruction lets continents meet and stop: the hand-built Artifexia project's
continental area never shrinks, and its colliding continents barely overlap. The core instead
resolves the overlap one sub-step's motion makes (up to some hundred km along the whole contact)
and keeps its volume as thicker crust, so continent is lost to every collision. With rigid
continents, before each sub-step, two continental plates converging across their boundary are
welded at once, where they are, when either continent is nearer that boundary than the sub-step
will close: they meet and from then on move as one.
"""
from __future__ import annotations

import numpy as np
import spherely

from paradise.core import sphere
from paradise.core.regions import CONTINENTAL_SHARE, SUTURE_KM
from paradise.core.topology import segment_motion

_KM_PER_S2_METRE = sphere.RADIUS_KM / spherely.EARTH_RADIUS_METERS
MARGIN_KM = 20.   # welded this far apart or less, whether or not they are converging fast
# Slower than this is no approach: near the pole of a rift's opening the motion is next to nothing,
# either way, and the pieces' continents still touch there.
MIN_CLOSING_CM_YR = .1
# Nor is a contact shorter than this: a kink in a new rift's boundary, beside the cut, where the
# pieces' continents still touch, welded them again in the step they rifted.
MIN_CONTACT_KM = 100.


def weld_on_approach(world, rotations: dict, dt: float) -> list:
    """Weld every pair of continental plates whose continents will touch within `dt`; returns events."""
    continents = {}

    def continent(plate):
        if plate not in continents:
            shapes = [c.shape for c in world.crust if c.plate == plate and c.kind == 'continent']
            continents[plate] = sphere.rotate_geography(sphere.union(*shapes), world.frame(plate), tolerance=1e-4) \
                if shapes else None
        return continents[plate]

    due = False
    for b in world.boundaries:
        if len(b.points) < 2 or not (world._continental(b.left) and world._continental(b.right)):
            continue
        if frozenset((b.left, b.right)) in getattr(world, 'kept_apart', ()):
            continue
        opening, _ = segment_motion(b.points, rotations.get(b.left, np.zeros(3)), rotations.get(b.right, np.zeros(3)))
        closing = np.maximum(-opening, 0.) * dt
        if not closing.any():
            continue
        a, c = continent(b.left), continent(b.right)
        if a is None or c is None:
            continue
        # Segment by segment: how near each continent is, and whether this segment closes that far
        # this sub-step. A whole boundary can close on its far side (round the back of the sphere,
        # from a rift's pole of opening) while the continents touch only where it opens. At a
        # trench it is the downgoing plate's continent that counts: it reaches the trench, and is
        # squeezed there, long before it reaches the other, while the overriding continent sits
        # near the trench by nature (an Andean margin), however wide the ocean still to close.
        lon, lat = sphere.to_lonlat(b.mids())
        probes = spherely.points(lon, lat)
        reach = closing + MARGIN_KM
        trench = [o for k, o in zip(b.kinds, b.overriding) if k == 'trench' and o in (b.left, b.right)]
        if trench and len(trench) >= len(b.kinds) / 2:
            over = max(set(trench), key=trench.count)
            down = c if over == b.left else a
            near = np.asarray(spherely.distance(probes, down)) * _KM_PER_S2_METRE <= reach
        else:
            near = (np.asarray(spherely.distance(probes, a)) * _KM_PER_S2_METRE <= reach) & \
                   (np.asarray(spherely.distance(probes, c)) * _KM_PER_S2_METRE <= reach)
        meeting = near & (closing > MIN_CLOSING_CM_YR * 10. * dt)   # 1 cm/yr is 10 km/Myr
        if float(np.sum(sphere.arc_km(b.points[:-1], b.points[1:])[meeting])) < MIN_CONTACT_KM:
            continue
        # Make the boundary a collision that has closed far enough to weld; the core does the rest.
        b.kinds = ['collision'] * len(b.kinds)
        b.overriding = [None] * len(b.kinds)
        b.closed_km = np.full(len(b.kinds), max(SUTURE_KM, getattr(world, 'suture_km', SUTURE_KM)))
        due = True
    if not due:
        return []
    events = world._sutures(rotations)
    for e in events:
        e['on_approach'] = True
    return events


__all__ = ['weld_on_approach', 'CONTINENTAL_SHARE']
