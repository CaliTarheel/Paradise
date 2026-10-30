"""Rigid continents: two continental plates weld where their continents meet, losing no continent."""
import unittest

import numpy as np

from paradise.core import sphere
from paradise.history import History, Project
from tests.test_paradise import EAST, four_plates, grid_lonlat


def two_continents():
    """Continents on plates 1 and 2, either side of their boundary along 0 degrees, 1,000 km apart."""
    plates, crust = four_plates()
    lon, lat = grid_lonlat()
    xyz = sphere.from_lonlat(lon, lat)
    for centre in ((-20., 0.), (20., 0.)):
        crust[sphere.arc_km(xyz, sphere.from_lonlat(*centre)[None, None]) < 2000.] = 1   # 12.6 Mkm2 each
    return plates, crust


def continent_km2(world):
    return sum(c.area_km2() for c in world.crust if c.kind == 'continent')


class RigidCollisionTests(unittest.TestCase):
    def run_collision(self, rigid):
        project = Project('test', *two_continents(), start_age_ma=1000., auto_trenches=True,
                          weld_continents_only=True, rigid_continents=rigid)
        history = History(project)
        before = continent_km2(history.current.world)
        welded = []
        for _ in range(4):
            snapshot = history.commit({1: EAST})
            welded += [e for entry in snapshot.entries for e in entry.get('events', []) if e['kind'] == 'suture']
        return before, continent_km2(history.current.world), welded, history.current.world

    def test_rigid_continents_weld_on_approach_and_keep_their_area(self):
        before, after, welded, world = self.run_collision(rigid=True)
        self.assertTrue(welded)
        self.assertTrue(all(e.get('on_approach') for e in welded))
        self.assertEqual(len({1, 2} & set(world.plates())), 1)
        self.assertGreater(after, .995 * before)   # nothing squeezed away

    def test_without_rigid_continents_the_collision_squeezes_continent(self):
        before, after, welded, world = self.run_collision(rigid=False)
        self.assertLess(after, .99 * before)


class ContinentalMarginTests(unittest.TestCase):
    def margin_closing(self, rigid):
        """Plate 1's continent reaches the boundary along 0 degrees; ocean plate 2 closes on it untrenched."""
        plates, crust = four_plates()
        lon, lat = grid_lonlat()
        crust[(lon > -30.) & (lon < 0.) & (np.abs(lat) < 40.)] = 1
        history = History(Project('test', plates, crust, start_age_ma=1000., auto_trenches=False,
                                  rigid_continents=rigid))
        before = continent_km2(history.current.world)
        history.commit({2: [90., 0., -.5]})
        return before, continent_km2(history.current.world)

    def test_rigid_continents_keep_a_margin_the_ocean_closes_on(self):
        before, after = self.margin_closing(rigid=True)
        self.assertGreater(after, .995 * before)

    def test_otherwise_the_margin_takes_half_the_overlap(self):
        before, after = self.margin_closing(rigid=False)
        self.assertLess(after, .99 * before)


if __name__ == '__main__':
    unittest.main()
