"""The Worldbuilding Pasta start and the advisor's proposals."""
import unittest

import numpy as np

from paradise import advisor as advice
from paradise.core import sphere
from paradise.history import History, Project, speed_cm_yr, pole_vector
from paradise.paint import world_from_grid
from paradise.starts import pasta_world


class PastaStartTests(unittest.TestCase):
    def test_one_ocean_plate_and_three_continental_plates(self):
        plates, crust = pasta_world()
        world, cratons = world_from_grid(plates, crust)
        self.assertEqual(world.plates(), [1, 2, 3, 4])
        self.assertAlmostEqual(sum(sphere.area_km2(r) for r in world.regions.values()) / sphere.SPHERE_KM2, 1., places=9)
        continent = {p: sum(c.area_km2() for c in world.crust if c.plate == p and c.kind == 'continent')
                     for p in world.plates()}
        self.assertEqual(continent[1], 0.)                       # the ocean plate
        for p in (2, 3, 4):
            self.assertGreater(continent[p], .5 * sphere.area_km2(world.regions[p]))   # continent and margin
        self.assertEqual(len(cratons), 3)


class AdvisorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        plates, crust = pasta_world()
        cls.history = History(Project('pasta', plates, crust, start_age_ma=1000., continents_override=True))

    def test_continental_plates_move_and_the_ocean_plate_barely_does(self):
        poles, events, reasons = advice.Advisor(3).propose(self.history)
        world = self.history.current.world
        for p, pole in poles.items():
            centre = advice._centre(world.regions[p])
            speed = speed_cm_yr(pole_vector(*pole), centre)
            if p == 1:
                self.assertLess(speed, 2.5)
            else:
                low, high = advice.CONTINENTAL_CM_YR
                self.assertGreater(speed, .5 * low)
                self.assertLess(speed, 1.01 * high)
        self.assertEqual(events, [])   # nothing converges yet, and nothing rifts in its first 120 Myr

    def test_trenches_put_the_continental_plate_over_the_ocean(self):
        a = advice.Advisor(3)
        poles, events, _ = a.propose(self.history)
        history = History(Project.from_json(self.history.project.to_json()))
        history.commit(poles, events=events)
        poles, events, _ = a.propose(history)
        subductions = [e for e in events if e['kind'] == 'subduction']
        history.commit(poles, events=events)
        for t in history.current.trenches:
            if 1 in t['plates']:
                self.assertNotEqual(t['overriding'], 1, 'the ocean plate overrides a continental one')
        for e in subductions:
            self.assertNotEqual(e['overriding'], 1)



class PangaeaTests(unittest.TestCase):
    def test_a_rifted_pangaea_comes_apart_and_stays_apart(self):
        from paradise.starts import pangaea_world
        plates, crust = pangaea_world()
        history = History(Project('pangaea', plates, crust, start_age_ma=1000., continents_override=True,
                                  weld_continents_only=True, rigid_continents=True))
        a = advice.Advisor(7, rift_at_start=True)
        a.rng = np.random.default_rng(0)
        advice_rift_chance, advice.RIFT_CHANCE = advice.RIFT_CHANCE, 1.
        try:
            poles, events, _ = a.propose(history)
        finally:
            advice.RIFT_CHANCE = advice_rift_chance
        self.assertEqual([e['kind'] for e in events if e['kind'] == 'rift'], ['rift'])
        snapshot = history.commit(poles, events=events)
        self.assertEqual(len(snapshot.world.plates()), 3)
        kinds = {k for b in snapshot.world.boundaries if {b.left, b.right} == {2, 3} for k in b.kinds}
        self.assertIn('ridge', kinds)
        self.assertNotIn('collision', kinds)


if __name__ == '__main__':
    unittest.main()
