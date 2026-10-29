"""Hotspots that leave trails, and large igneous provinces."""
import math
import os
import tempfile
import unittest

import numpy as np
import pygplates

from paradise import gplates, plumes
from paradise.core import sphere
from paradise.history import History, Project
from tests.test_paradise import EAST, four_plates


def continent_km2(world, plate):
    return sum(c.area_km2() for c in world.crust if c.kind == 'continent' and c.plate == plate)


class HotspotTests(unittest.TestCase):
    def setUp(self):
        self.project = Project('test', *four_plates(), start_age_ma=1000., auto_trenches=False)

    def test_a_hotspot_leaves_a_trail_that_rides_away_with_its_plate(self):
        history = History(self.project)
        history.commit({1: EAST}, events=[dict(kind='hotspot', point=[-60., 0.], name='Hawaii')])
        snapshot = history.commit({1: EAST})
        hotspot = snapshot.hotspots[0]
        self.assertEqual(hotspot['name'], 'Hawaii')
        self.assertEqual({m['plate'] for m in hotspot['marks']}, {1})
        self.assertEqual(len(hotspot['marks']), 8)   # one per 2.5 Myr sub-step over two steps
        world = snapshot.world
        # The oldest mark has moved east with plate 1 by the rotation since it was made: 0.5 degrees/Myr
        # about the north pole for 17.5 Myr.
        oldest = sphere.rotate(np.array(hotspot['marks'][0]['point']), world.frame(1))
        lon, lat = sphere.to_lonlat(oldest)
        self.assertAlmostEqual(float(lon), -60. + .5 * 17.5, delta=.01)
        self.assertAlmostEqual(float(lat), 0., delta=.01)
        trail = plumes.trails(snapshot.hotspots)
        self.assertEqual(len(trail), 1)

    def test_an_ended_hotspot_leaves_no_more_marks(self):
        history = History(self.project)
        history.commit({1: EAST}, events=[dict(kind='hotspot', point=[-60., 0.])])
        history.commit({1: EAST}, events=[dict(kind='end_hotspot', hotspot=1)])
        self.assertEqual(len(history.current.hotspots[0]['marks']), 4)
        self.assertEqual(history.current.hotspots[0]['ended_myr'], 10.)

    def test_hotspots_and_trails_export(self):
        history = History(self.project)
        history.commit({1: EAST}, events=[dict(kind='hotspot', point=[-60., 0.], name='Hawaii')])
        history.commit({1: EAST})
        with tempfile.TemporaryDirectory() as directory:
            gplates.export(history, directory)
            spots = list(pygplates.FeatureCollection(os.path.join(directory, 'Hotspots.gpml')))
            self.assertEqual(len(spots), 1)
            self.assertEqual(spots[0].get_reconstruction_plate_id(), 0)
            trails = list(pygplates.FeatureCollection(os.path.join(directory, 'Hotspot trails.gpml')))
            self.assertGreater(len(trails), 0)
            self.assertTrue(all(f.get_reconstruction_plate_id() == 1 for f in trails))


class LipTests(unittest.TestCase):
    def test_a_lip_over_ocean_makes_a_plateau_and_leaves_a_hotspot(self):
        project = Project('test', *four_plates(), start_age_ma=1000., auto_trenches=False)
        history = History(project)
        snapshot = history.commit({}, events=[dict(kind='lip', point=[-120., 20.], radius_km=600., name='Ontong')])
        plateau = continent_km2(snapshot.world, 1)
        self.assertAlmostEqual(plateau, math.pi * 600. ** 2, delta=.03 * math.pi * 600. ** 2)
        self.assertEqual([c.thickness_km for c in snapshot.world.crust if c.kind == 'continent'], [plumes.PLATEAU_KM])
        lips = [f for f in snapshot.features if f.kind == 'lip']
        self.assertEqual(len(lips), 1)
        self.assertEqual([h['name'] for h in snapshot.hotspots], ['Ontong'])
        events = {e['kind'] for e in snapshot.entries[0]['events']}
        self.assertIn('lip_erupted', events)
        with tempfile.TemporaryDirectory() as directory:
            gplates.export(history, directory)
            active = list(pygplates.FeatureCollection(os.path.join(directory, 'Active LIPs.gpml')))
            self.assertEqual(len(active), 1)
            self.assertEqual(active[0].get_feature_type(), pygplates.FeatureType.gpml_large_igneous_province)

    def test_a_lip_without_a_hotspot(self):
        history = History(Project('test', *four_plates(), start_age_ma=1000.))
        snapshot = history.commit({}, events=[dict(kind='lip', point=[-120., 20.], with_hotspot=False)])
        self.assertEqual(snapshot.hotspots, [])


if __name__ == '__main__':
    unittest.main()
