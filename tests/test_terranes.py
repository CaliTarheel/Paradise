"""Island arcs behind trenches, and terranes accreting instead of being subducted."""
import os
import tempfile
import unittest

import numpy as np
import pygplates

from paradise import gplates, terranes
from paradise.core import sphere
from paradise.history import History, Project
from tests.test_paradise import EAST, STRETCH, centroid, four_plates, grid_lonlat


def world_with_a_microcontinent(radius_km=350.):
    """four_plates, with a small continent on plate 1 15 degrees west of the trench line at 0."""
    plates, crust = four_plates()
    lon, lat = grid_lonlat()
    distance = sphere.arc_km(sphere.from_lonlat(lon, lat), sphere.from_lonlat(-15., 5.)[None, None])
    crust[distance < radius_km] = 1
    return Project('test', plates, crust, start_age_ma=1000., auto_trenches=False)


TRENCH = dict(kind='subduction', plates=[1, 2], overriding=2, points=STRETCH)   # plate 1 goes under


def continent_km2(world, plate):
    return sum(c.area_km2() for c in world.crust if c.kind == 'continent' and c.plate == plate)


class ArcTests(unittest.TestCase):
    def test_an_arc_rises_behind_a_trench_that_has_lasted_a_step(self):
        history = History(world_with_a_microcontinent())
        first = history.commit({1: EAST}, events=[TRENCH])
        self.assertEqual(first.features, [])   # a trench declared this step has fed no arc yet
        second = history.commit({1: EAST})
        arcs = [f for f in second.features if f.kind == 'island_arc']
        self.assertEqual(len(arcs), 1)
        arc = arcs[0]
        self.assertEqual(arc.plate, 2)
        self.assertAlmostEqual(continent_km2(second.world, 2), arc.area_km2(), delta=1.)
        # The arc lies on the overriding plate, ARC_OFFSET_KM inboard of the trench (east of 0 degrees).
        lon, lat = sphere.to_lonlat(sphere.rotate(centroid(arc.shape), second.world.frame(2)))
        trench_lon = np.mean([p[0] for t in second.trenches for p in t['points'] if abs(p[1] - float(lat)) < 20.])
        offset_km = (float(lon) - trench_lon) * np.pi / 180. * sphere.RADIUS_KM * np.cos(np.radians(float(lat)))
        self.assertAlmostEqual(offset_km, terranes.ARC_OFFSET_KM, delta=60.)
        self.assertIn('arc_born', {e['kind'] for e in second.entries[-1]['events']})

    def test_the_offset_line_lies_on_the_side_asked(self):
        line = sphere.from_lonlat(np.zeros(5), np.linspace(-20., 20., 5))   # northwards along 0 degrees
        left = terranes.offset_line(line, 100., to_left=True)
        right = terranes.offset_line(line, 100., to_left=False)
        self.assertTrue(np.all(sphere.to_lonlat(left)[0] < 0.))    # west is on the left going north
        self.assertTrue(np.all(sphere.to_lonlat(right)[0] > 0.))
        np.testing.assert_allclose(sphere.arc_km(left, line), 100., atol=.5)


class TerraneTests(unittest.TestCase):
    def test_a_microcontinent_accretes_and_the_trench_steps_outboard(self):
        history = History(world_with_a_microcontinent())
        history.commit({1: EAST}, events=[TRENCH])
        history.commit({1: EAST})
        before = continent_km2(history.current.world, 1)
        snapshot = history.commit({1: EAST})
        events = [e for entry in snapshot.entries for e in entry.get('events', []) if e['kind'] == 'accreted']
        self.assertEqual(len(events), 1)
        self.assertEqual((events[0]['plate'], events[0]['into']), (1, 2))
        self.assertAlmostEqual(events[0]['area_km2'], before, delta=1.)
        world = snapshot.world
        self.assertEqual(continent_km2(world, 1), 0.)
        terrane = [f for f in snapshot.features if f.kind == 'terrane']
        self.assertEqual(len(terrane), 1)
        self.assertEqual(terrane[0].plate, 2)
        # It was neither subducted nor collided: its area is all there, and the boundary stayed a trench.
        self.assertAlmostEqual(terrane[0].area_km2(), before, delta=.01 * before)
        # Along the declared stretch near 0 degrees (its ends, beyond 50 degrees, were never declared;
        # the other boundary between the two plates, a ridge, is at the date line).
        kinds = set()
        for b in world.boundaries:
            if (b.left, b.right) == (1, 2):
                lon, lat = sphere.to_lonlat(b.mids())
                kinds |= {k for k, x, y in zip(b.kinds, lon, lat) if abs(y) < 45. and abs(x) < 30.}
        self.assertEqual(kinds, {'trench'})
        # The trench now runs round the outboard side of the terrane (west of it).
        west = min(p[0] for t in snapshot.trenches for p in t['points'] if abs(p[1] - 5.) < 3.)
        lon, _ = sphere.to_lonlat(sphere.rotate(centroid(terrane[0].shape), world.frame(2)))
        self.assertLess(west, float(lon))
        self.assertEqual(terrane[0].accreted[0]['to_plate'], 2)

    def test_a_continent_larger_than_a_terrane_does_not_accrete(self):
        project = world_with_a_microcontinent(radius_km=900.)   # about 2.5 million km2
        self.assertGreater(continent_km2(History(project).current.world, 1), terranes.TERRANE_MAX_KM2)
        history = History(project)
        history.commit({1: EAST}, events=[TRENCH])
        history.commit({1: EAST})
        snapshot = history.commit({1: EAST})
        kinds = {e['kind'] for entry in snapshot.entries for e in entry.get('events', [])}
        self.assertNotIn('accreted', kinds)

    def test_arcs_and_terranes_export_to_their_layers(self):
        project = world_with_a_microcontinent()
        history = History(project)
        history.commit({1: EAST}, events=[TRENCH])
        history.commit({1: EAST})
        history.commit({1: EAST})
        with tempfile.TemporaryDirectory() as directory:
            gplates.export(history, directory)
            arcs = pygplates.FeatureCollection(os.path.join(directory, 'Island arcs.gpml'))
            accreted = pygplates.FeatureCollection(os.path.join(directory, 'Accreted terranes.gpml'))
            self.assertGreater(len(arcs), 0)
            self.assertGreater(len(accreted), 0)
            self.assertTrue(all(f.get_feature_type() == pygplates.FeatureType.gpml_island_arc for f in arcs))
            description = list(accreted)[0].get_description()
            self.assertIn('accreted from plate 1 to plate 2', description)


if __name__ == '__main__':
    unittest.main()
