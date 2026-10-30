"""Declared rifts and back-arcs: a plate split in two at the start of a step."""
import os
import tempfile
import unittest

import numpy as np
import pygplates

from paradise import gplates
from paradise.core import sphere
from paradise.history import History, Project
from tests.test_paradise import STRETCH, four_plates, grid_lonlat

# Across plate 1 of four_plates (west of 0, between 60 S and 60 N), along 90 W through its continent.
LINE = [[-90., -70.], [-90., 70.]]
AWAY = [90., 0., .4]   # the new plate, east of the line, turns east about the north pole


def with_a_continent():
    plates, crust = four_plates()
    lon, lat = grid_lonlat()
    distance = sphere.arc_km(sphere.from_lonlat(lon, lat), sphere.from_lonlat(-90., 0.)[None, None])
    crust[distance < 2500.] = 1
    return plates, crust


def continent_km2(world, plate):
    return sum(c.area_km2() for c in world.crust if c.kind == 'continent' and c.plate == plate)


class RiftTests(unittest.TestCase):
    def setUp(self):
        self.project = Project('test', *with_a_continent(), start_age_ma=1000., auto_trenches=False)

    def test_a_rift_splits_the_plate_and_its_continent_and_opens_a_ridge(self):
        history = History(self.project)
        continent = continent_km2(history.current.world, 1)
        snapshot = history.commit({}, events=[dict(kind='rift', plate=1, points=LINE, new_side='right', pole=AWAY)])
        world = snapshot.world
        self.assertEqual(world.plates(), [1, 2, 3, 4, 5])
        self.assertEqual(world.parents[5], 1)
        rifts = [e for e in snapshot.entries[0]['events'] if e['kind'] == 'rift']
        self.assertEqual([(e['plate'], e['new']) for e in rifts], [(1, 5)])
        # Each side keeps its share of the continent, and none is lost.
        self.assertGreater(continent_km2(world, 1), .3 * continent)
        self.assertGreater(continent_km2(world, 5), .3 * continent)
        self.assertAlmostEqual(continent_km2(world, 1) + continent_km2(world, 5), continent, delta=.01 * continent)
        kinds = {k for b in world.boundaries if (b.left, b.right) == (1, 5) for k in b.kinds}
        self.assertIn('ridge', kinds)
        self.assertGreater(sum(e['gained']['ridge'] for e in snapshot.entries), 1e5)
        self.assertEqual([(r['plate'], r['new']) for r in snapshot.rifts], [(1, 5)])
        # The new plate's pole is kept in the log, and the history replays the same.
        replayed = History(Project.from_json(self.project.to_json()))
        self.assertEqual(replayed.current.world.plates(), [1, 2, 3, 4, 5])
        self.assertAlmostEqual(continent_km2(replayed.current.world, 5), continent_km2(world, 5), delta=1.)

    def test_a_rift_without_its_own_pole_moves_with_its_parent(self):
        history = History(self.project)
        history.commit({1: [0., 0., .3]}, events=[dict(kind='rift', plate=1, points=LINE, new_side='right')])
        frames = history.current.world.frames
        np.testing.assert_allclose(frames[5], frames[1], atol=1e-12)

    def test_a_line_that_does_not_cross_the_plate_is_refused(self):
        history = History(self.project)
        snapshot = history.preview({}, [dict(kind='rift', plate=1, points=[[-100., 10.], [-80., 10.]], new_side='left')])
        self.assertIn('split_refused', {w['kind'] for w in snapshot.warnings})
        self.assertEqual(snapshot.world.plates(), [1, 2, 3, 4])

    def test_a_rifted_plate_rotates_relative_to_its_parent_in_the_export(self):
        history = History(self.project)
        history.commit({1: [0., 0., .3]}, events=[dict(kind='rift', plate=1, points=LINE, new_side='right', pole=AWAY)])
        history.commit({1: [0., 0., .3], 5: AWAY})
        rows = [r for r in gplates.rotation_rows(history) if r[0] == 5]
        self.assertTrue(rows)
        self.assertEqual({r[5] for r in rows}, {1})
        with tempfile.TemporaryDirectory() as directory:
            gplates.export(history, directory)
            model = pygplates.RotationModel(os.path.join(directory, 'rotation.rot'))
            world = history.current.world
            age = self.project.age_ma(history.current.time_myr)
            # The same total rotation, through the tree.
            point = sphere.from_lonlat(-60., 10.)
            moved = model.get_rotation(age, 5) * pygplates.PointOnSphere(10., -60.)
            np.testing.assert_allclose(np.array(moved.to_xyz()), sphere.rotate(point, world.frame(5)), atol=1e-6)

    def test_the_rift_line_exports(self):
        history = History(self.project)
        history.commit({}, events=[dict(kind='rift', plate=1, points=LINE, new_side='right', pole=AWAY)])
        with tempfile.TemporaryDirectory() as directory:
            gplates.export(history, directory)
            rifts = pygplates.FeatureCollection(os.path.join(directory, 'Rifts.gpml'))
            self.assertEqual(len(rifts), 1)
            self.assertEqual(list(rifts)[0].get_feature_type(), pygplates.FeatureType.gpml_continental_rift)

    def test_a_craton_goes_with_its_plate_when_the_plate_is_welded_away(self):
        plates, crust = with_a_continent()
        lon, lat = grid_lonlat()
        # A craton on the northern cap, plate 3, which is smaller than plate 1 and so is the one welded away.
        crust[sphere.arc_km(sphere.from_lonlat(lon, lat), sphere.from_lonlat(-90., 75.)[None, None]) < 500.] = 2
        history = History(Project('test', plates, crust, start_age_ma=1000.))
        self.assertEqual([p for p, _ in history.current.cratons], [3])
        snapshot = history.commit({}, events=[dict(kind='weld', plates=[1, 3])])
        self.assertNotIn(3, snapshot.world.plates())
        self.assertEqual([p for p, _ in snapshot.cratons], [1])

    def test_a_rift_welded_back_exports_as_failed(self):
        history = History(self.project)
        history.commit({}, events=[dict(kind='rift', plate=1, points=LINE, new_side='right', pole=AWAY)])
        history.commit({}, events=[dict(kind='weld', plates=[1, 5])])
        with tempfile.TemporaryDirectory() as directory:
            gplates.export(history, directory)
            failed = pygplates.FeatureCollection(os.path.join(directory, 'Failed rifts.gpml'))
            active = pygplates.FeatureCollection(os.path.join(directory, 'Rifts.gpml'))
            self.assertEqual(len(failed), 1)   # the last step, after the weld
            self.assertEqual(len(active), 1)   # the step it opened

    def test_old_ocean_detaches_from_its_plate(self):
        project = Project('test', *with_a_continent(), start_age_ma=1000., ocean_age_myr=200.)
        history = History(project)
        snapshot = history.commit({}, events=[dict(kind='detach', plate=1, min_age_myr=150., pole=[0., 0., .3])])
        world = snapshot.world
        self.assertEqual(world.plates(), [1, 2, 3, 4, 5])
        self.assertEqual(world.parents[5], 1)
        # All the old ocean went; the continent and its plate stayed.
        self.assertEqual(sum(c.area_km2() for c in world.crust if c.plate == 5 and c.kind == 'continent'), 0.)
        self.assertGreater(continent_km2(world, 1), 0.)
        self.assertGreater(sphere.area_km2(world.regions[5]), 1e7)


class BackArcTests(unittest.TestCase):
    def test_a_backarc_rolls_back_and_opens_ocean_on_the_overriding_side(self):
        project = Project('test', *four_plates(), start_age_ma=1000., auto_trenches=False)
        history = History(project)
        # Plate 2 overrides plate 1 along 0 degrees; plate 1 stands still.
        first = history.commit({}, events=[dict(kind='subduction', plates=[1, 2], overriding=2, points=STRETCH)])
        trench = first.trenches[0]
        # The forearc rolls back west, towards the trench, about the north pole.
        snapshot = history.commit({}, events=[dict(kind='backarc', plates=[1, 2], overriding=2,
                                                   points=trench['points'], pole=[90., 0., -.4])])
        world = snapshot.world
        self.assertEqual(world.plates(), [1, 2, 3, 4, 5])
        forearc = 5
        self.assertEqual(world.parents[forearc], 2)
        # The forearc now overrides the trench, and a ridge opens between it and plate 2 behind it.
        self.assertEqual({t['overriding'] for t in snapshot.trenches}, {forearc})
        kinds = {k for b in world.boundaries if (b.left, b.right) == (2, forearc) for k in b.kinds}
        self.assertIn('ridge', kinds)
        self.assertGreater(sum(e['gained']['ridge'] for e in snapshot.entries), 1e5)
        # About BACKARC_KM wide all along the trench.
        area = sphere.area_km2(world.regions[forearc])
        self.assertGreater(area, .5 * 250. * trench['length_km'])


if __name__ == '__main__':
    unittest.main()
