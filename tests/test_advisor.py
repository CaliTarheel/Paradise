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


def cap(lon, lat, radius_km):
    centre = sphere.from_lonlat(lon, lat)
    east = sphere.unit(np.cross([0., 0., 1.], centre))
    north = np.cross(centre, east)
    a = radius_km / sphere.RADIUS_KM
    t = np.linspace(0., 2. * np.pi, 72, endpoint=False)
    return sphere.polygon(np.cos(a) * centre + np.sin(a) * (np.outer(np.cos(t), east) + np.outer(np.sin(t), north)))


class OpeningAxisTests(unittest.TestCase):
    def test_a_rift_opens_wherever_it_crosses_the_continent(self):
        import math
        # Two continents on the equator's great circle, 100 degrees apart: the rift along it must
        # open across both, so its pole goes in the 260 degrees of ocean.
        land = sphere.union(*[cap(lon, 0., 1500.) for lon in (0., 100.)])
        centre, heading = sphere.from_lonlat(0., 0.), sphere.from_lonlat(90., 0.)
        axis = advice._opening_axis(land, centre, heading)
        side = np.cross(centre, heading)
        for lon in np.arange(-13., 114., 1.):
            point = sphere.from_lonlat(lon, 0.)
            self.assertGreater(np.dot(np.cross(axis, point), side), 0., lon)
        self.assertAlmostEqual(float(np.dot(axis, side)), 0.)

    def test_none_when_the_continent_spans_more_than_half_the_circle(self):
        land = sphere.union(*[cap(lon, 0., 1500.) for lon in (0., 100., -100.)])
        self.assertIsNone(advice._opening_axis(land, sphere.from_lonlat(0., 0.), sphere.from_lonlat(90., 0.)))


class ClockTests(unittest.TestCase):
    def snapshot(self, time_myr, *events):
        class Snapshot:
            pass
        s = Snapshot()
        s.time_myr, s.entries = time_myr, [dict(events=list(events))]
        return s

    def test_a_detached_old_ocean_does_not_reset_the_rift_clock(self):
        a = advice.Advisor(1)
        a.committed(self.snapshot(50., dict(kind='rift', plate=3, new=9, name='Plate 3 old ocean')))
        self.assertNotIn(3, a.last_rift)
        self.assertEqual(a.last_rift[9], 50.)
        a.pending = {3: np.zeros(3)}
        a.committed(self.snapshot(60., dict(kind='rift', plate=3, new=10, name='Plate 3 rift')))
        self.assertEqual(a.last_rift[3], 60.)
        self.assertIn(10, a.spreading)

    def test_the_clocks_are_read_back_from_a_history(self):
        class History:
            pass
        h = History()
        h.project = History()
        h.project.steps = [dict(poles={'1': [0., 0., .2]}, events=[dict(kind='rift', plate=1, pole=[0., 0., .5])]),
                           dict(poles={}, events=[dict(kind='detach', plate=2)])]
        h.snapshots = [None, self.snapshot(10., dict(kind='rift', plate=1, new=2, name='Plate 1 rift')),
                       self.snapshot(20., dict(kind='rift', plate=2, new=3, name='Plate 2 old ocean'))]
        a = advice.Advisor(1)
        a.remember(h)
        self.assertEqual(a.last_rift, {1: 10., 2: 10., 3: 20.})
        self.assertEqual(a.last_detach, {2: 20.})
        parent, away, until = a.spreading[2]
        self.assertEqual((parent, until), (1, 10. + advice.SEPARATION_MYR))
        np.testing.assert_allclose(away, pole_vector(0., 0., .3), atol=1e-12)


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
