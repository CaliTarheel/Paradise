"""Painted worlds, the step history and the GPlates export."""
import os
import tempfile
import unittest

import numpy as np
import pygplates
import spherely

from paradise import gplates, paint
from paradise.core import sphere
from paradise.history import History, Project, pole_vector, vector_pole
from paradise.paint import CONTINENT, CRATON, cell_km, world_from_grid

H, W = 90, 180


def grid_lonlat():
    lat = 90. - (np.arange(H) + .5) * 180. / H
    lon = -180. + (np.arange(W) + .5) * 360. / W
    return np.meshgrid(lon, lat)


def three_plates():
    """Plate 1 west, 2 east, 3 a polar cap; a continent (with a craton) on plate 1 at 90 W."""
    lon, lat = grid_lonlat()
    plates = np.where(lon < 0., 1, 2)
    plates[lat > 60.] = 3
    distance = sphere.arc_km(sphere.from_lonlat(lon, lat), sphere.from_lonlat(-90., 0.)[None, None])
    crust = np.where(distance < 2500., CONTINENT, 0)
    crust[distance < 1000.] = CRATON
    return plates, crust


def painted_area_km2(mask):
    lat = np.radians(90. - (np.arange(H + 1)) * 180. / H)
    band = (np.sin(lat[:-1]) - np.sin(lat[1:])) * 2. * np.pi / W * sphere.RADIUS_KM ** 2
    return float((mask * band[:, None]).sum())


def centroid(geography):
    c = spherely.centroid(geography)
    return sphere.from_lonlat(spherely.get_x(c), spherely.get_y(c))


class PaintTests(unittest.TestCase):
    def test_plates_tile_the_sphere_and_keep_their_painted_areas(self):
        plates, crust = three_plates()
        world, cratons = world_from_grid(plates, crust)
        self.assertEqual(world.plates(), [1, 2, 3])
        total = sum(sphere.area_km2(r) for r in world.regions.values())
        self.assertAlmostEqual(total / sphere.SPHERE_KM2, 1., places=6)
        # Taking the staircase out moves area between neighbours by about a cell per corner, no more.
        tolerance = 3. * cell_km(W) ** 2 * 40
        for p in (1, 2, 3):
            self.assertLess(abs(sphere.area_km2(world.regions[p]) - painted_area_km2(plates == p)), tolerance)
        continent = sum(c.area_km2() for c in world.crust if c.kind == 'continent')
        self.assertLess(abs(continent - painted_area_km2(crust > 0)), tolerance)
        self.assertEqual(len(cratons), 1)
        self.assertEqual(cratons[0][0], 1)

    def test_a_plate_larger_than_a_hemisphere_and_one_across_the_date_line(self):
        lon, lat = grid_lonlat()
        plates = np.ones((H, W), np.int64)
        plates[(np.abs(lon) > 150.) & (np.abs(lat) < 30.)] = 2   # straddles 180 degrees
        world, _ = world_from_grid(plates, np.zeros((H, W), np.int64))
        self.assertEqual(world.plates(), [1, 2])
        self.assertGreater(sphere.area_km2(world.regions[1]), .5 * sphere.SPHERE_KM2)
        lon_c, lat_c = sphere.to_lonlat(centroid(world.regions[2]))
        self.assertGreater(abs(float(lon_c)), 170.)
        self.assertLess(abs(float(lat_c)), 1.)

    def test_coastline_detail_keeps_the_tiling_and_changes_the_coast(self):
        plates, crust = three_plates()
        plain, _ = world_from_grid(plates, crust)
        detailed, _ = world_from_grid(plates, crust, detail_km=paint.DETAIL_KM)
        total = sum(sphere.area_km2(r) for r in detailed.regions.values())
        self.assertAlmostEqual(total / sphere.SPHERE_KM2, 1., places=9)
        for p, region in detailed.regions.items():
            crust_km2 = sum(c.area_km2() for c in detailed.crust if c.plate == p)
            self.assertAlmostEqual(crust_km2, sphere.area_km2(region), delta=1.)
        coast = [c for c in detailed.crust if c.kind == 'continent'][0]
        before = [c for c in plain.crust if c.kind == 'continent'][0]
        moved = sphere.area_km2(sphere.difference(coast.shape, before.shape)) + \
            sphere.area_km2(sphere.difference(before.shape, coast.shape))
        self.assertGreater(moved, .02 * before.area_km2())   # the coast really moved
        self.assertLess(abs(coast.area_km2() - before.area_km2()), .15 * before.area_km2())

    def test_plates_meeting_at_a_pole(self):
        lon, lat = grid_lonlat()
        plates = 1 + ((lon + 180.) // 120.).astype(np.int64)   # three wedges from pole to pole
        world, _ = world_from_grid(plates, np.zeros((H, W), np.int64))
        self.assertEqual(world.plates(), [1, 2, 3])
        for p in (1, 2, 3):
            self.assertAlmostEqual(sphere.area_km2(world.regions[p]) / sphere.SPHERE_KM2, 1. / 3., places=3)


class HistoryTests(unittest.TestCase):
    def setUp(self):
        plates, crust = three_plates()
        self.project = Project('test', plates, crust, start_age_ma=1000.)

    def test_poles_round_trip(self):
        lat, lon, rate = vector_pole(pole_vector(-20., 140., .7))
        self.assertAlmostEqual(lat, -20.)
        self.assertAlmostEqual(lon, 140.)
        self.assertAlmostEqual(rate, .7)

    def test_a_step_moves_a_plate_and_opens_a_ridge_behind_it(self):
        history = History(self.project)
        before = centroid([c for c in history.current.world.crust if c.kind == 'continent'][0].shape)
        snapshot = history.commit({1: [0., 0., .5]}, 'east')
        self.assertAlmostEqual(snapshot.time_myr, 10.)
        world = snapshot.world
        piece = [c for c in world.crust if c.kind == 'continent'][0]
        after = sphere.rotate(centroid(piece.shape), world.frame(piece.plate))
        # 5 degrees about a pole at (0, 0): the continent at 90 W on the equator goes 5 degrees south.
        self.assertAlmostEqual(float(sphere.to_lonlat(after)[1]), sphere.to_lonlat(before)[1] - 5., delta=.2)
        kinds = {k for b in world.boundaries for k in b.kinds}
        self.assertIn('ridge', kinds)
        self.assertIn('trench', kinds)
        self.assertEqual(self.project.steps[0]['note'], 'east')
        for entry in snapshot.entries:
            self.assertLess(abs(entry['tiling_error_km2']), 1e-4 * sphere.SPHERE_KM2)

    def test_going_back_and_replaying_give_the_same_world(self):
        history = History(self.project)
        history.commit({1: [0., 0., .5]})
        history.commit({1: [0., 0., .5], 2: [30., 100., -.3]})
        areas = {p: sphere.area_km2(r) for p, r in history.current.world.regions.items()}
        replayed = History(Project.from_json(self.project.to_json()))
        again = {p: sphere.area_km2(r) for p, r in replayed.current.world.regions.items()}
        self.assertEqual(sorted(areas), sorted(again))
        for p in areas:
            self.assertAlmostEqual(areas[p], again[p], delta=1.)
        history.back_to(1)
        self.assertEqual(len(history.snapshots), 2)
        self.assertEqual(len(self.project.steps), 1)

    def test_the_cache_is_used_only_for_the_same_steps(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = os.path.join(directory, 'world.cache')
            history = History(self.project, cache)
            history.commit({1: [0., 0., .5]})
            history.save_cache()
            reloaded = History(Project.from_json(self.project.to_json()), cache)
            self.assertIs(type(reloaded.snapshots[1].world), type(history.snapshots[1].world))
            self.assertAlmostEqual(reloaded.current.time_myr, 10.)
            changed = Project.from_json(self.project.to_json())
            changed.steps[0]['poles'] = {'1': [0., 0., .25]}
            other = History(changed, cache)
            piece = [c for c in other.current.world.crust if c.kind == 'continent'][0]
            moved = sphere.rotate(centroid(piece.shape), other.current.world.frame(piece.plate))
            self.assertAlmostEqual(float(sphere.to_lonlat(moved)[1]), -2.5, delta=.2)


def four_plates():
    """Plates 1 (west) and 2 (east) meeting along 0 and 180 degrees, with polar caps 3 and 4."""
    lon, lat = grid_lonlat()
    plates = np.where(lon < 0., 1, 2)
    plates[lat > 60.] = 3
    plates[lat < -60.] = 4
    return plates, np.zeros((H, W), np.int64)


EAST = [90., 0., .5]           # plate 1 turns east about the north pole: closing on plate 2 along 0 degrees
STRETCH = [[0., -50.], [0., 50.]]


class DeclaredSubductionTests(unittest.TestCase):
    def setUp(self):
        self.project = Project('test', *four_plates(), start_age_ma=1000., auto_trenches=False)

    def test_without_a_declaration_no_trench_forms(self):
        snapshot = History(self.project).commit({1: EAST})
        self.assertEqual(snapshot.trenches, [])
        self.assertNotIn('trench', {k for b in snapshot.world.boundaries for k in b.kinds})
        self.assertGreater(sum(e['lost']['compression'] for e in snapshot.entries), 1e6)

    def test_a_declared_trench_consumes_the_plate_that_goes_under(self):
        history = History(self.project)
        snapshot = history.commit({1: EAST}, events=[dict(kind='subduction', plates=[2, 1], overriding=1,
                                                         points=STRETCH)])
        self.assertEqual(len(snapshot.trenches), 1)
        trench = snapshot.trenches[0]
        self.assertEqual((trench['overriding'], trench['downgoing']), (1, 2))
        # 0.5 degrees/Myr about the pole is 5.6 cm/yr at the equator, less towards the junctions at 60 degrees.
        self.assertGreater(trench['closing_cm_yr'], 3.)
        self.assertLess(trench['closing_cm_yr'], 5.6)
        self.assertAlmostEqual(trench['slab_age_myr'], 70., delta=1.)
        subducted = sum(e['lost']['subducted'] for e in snapshot.entries)
        self.assertGreater(subducted, .8 * trench['consumed_km2_per_step'])
        self.assertLess(subducted, 1.2 * trench['consumed_km2_per_step'])
        self.assertEqual(self.project.steps[0]['events'][0]['overriding'], 1)
        # The trench lasts into the next step without being declared again, and replays the same.
        again = history.commit({1: EAST})
        self.assertEqual([(t['overriding'], t['downgoing']) for t in again.trenches], [(1, 2)])
        replayed = History(Project.from_json(self.project.to_json()))
        self.assertEqual([t['length_km'] for t in replayed.current.trenches], [t['length_km'] for t in again.trenches])

    def test_declaring_again_flips_the_polarity_and_ending_stops_it(self):
        history = History(self.project)
        first = history.commit({1: EAST}, events=[dict(kind='subduction', plates=[1, 2], overriding=1, points=STRETCH)])
        flipped = history.commit({1: EAST}, events=[dict(kind='subduction', plates=[1, 2], overriding=2,
                                                        points=first.trenches[0]['points'])])
        self.assertEqual([(t['overriding'], t['downgoing']) for t in flipped.trenches], [(2, 1)])
        ended = history.commit({1: EAST}, events=[dict(kind='end_subduction', plates=[1, 2],
                                                      points=flipped.trenches[0]['points'])])
        self.assertEqual(ended.trenches, [])

    def test_a_trench_whose_plates_pull_apart_is_reported(self):
        history = History(self.project)
        west = [90., 0., -.5]
        snapshot = history.commit({1: west}, events=[dict(kind='subduction', plates=[1, 2], overriding=1, points=STRETCH)])
        self.assertLess(snapshot.trenches[0]['closing_cm_yr'], 0.)
        self.assertIn('trench_opening', {w['kind'] for w in snapshot.warnings})
        self.assertGreater(sum(e['gained']['trench_opening'] for e in snapshot.entries), 1e6)

    def test_bad_declarations_are_refused(self):
        history = History(self.project)
        with self.assertRaises(ValueError):
            history.preview({}, [dict(kind='subduction', plates=[1, 2], overriding=3, points=STRETCH)])
        with self.assertRaises(ValueError):
            history.preview({}, [dict(kind='subduction', plates=[1, 1], overriding=1, points=STRETCH)])
        missed = history.preview({}, [dict(kind='subduction', plates=[3, 4], overriding=3, points=STRETCH)])
        self.assertIn('declaration_missed', {w['kind'] for w in missed.warnings})


class ExportTests(unittest.TestCase):
    def test_gplates_reconstructs_what_the_model_shows(self):
        plates, crust = three_plates()
        project = Project('test', plates, crust, start_age_ma=1000.)
        history = History(project)
        history.commit({1: [0., 0., .5]})
        history.commit({1: [0., 0., .5], 2: [30., 100., -.3]})
        with tempfile.TemporaryDirectory() as directory:
            files = gplates.export(history, directory)
            names = {os.path.basename(f) for f in files}
            self.assertTrue({'rotation.rot', 'Plates.gpml', 'Continents.gpml', 'Cratons.gpml', 'Ocean crust.gpml',
                             'Divergent.gpml', 'Convergent.gpml', 'Transform.gpml'} <= names)
            model = pygplates.RotationModel(os.path.join(directory, 'rotation.rot'))
            for snapshot in history.snapshots:
                age = project.age_ma(snapshot.time_myr)
                reconstructed = []
                pygplates.reconstruct(os.path.join(directory, 'Continents.gpml'), model, reconstructed, age)
                # The painted continent; island arcs behind the trenches are continental crust too.
                world = snapshot.world
                piece = max((c for c in world.crust if c.kind == 'continent'), key=lambda c: c.area_km2())
                got = max((r.get_reconstructed_geometry() for r in reconstructed), key=lambda g: g.get_area())
                want = sphere.rotate(centroid(piece.shape), world.frame(piece.plate))
                lat, lon = got.get_interior_centroid().to_lat_lon()
                self.assertLess(float(sphere.arc_km(sphere.from_lonlat(lon, lat), want)), 5.)
                self.assertAlmostEqual(got.get_area() * sphere.RADIUS_KM ** 2 / piece.area_km2(), 1., places=3)
            divergent = pygplates.FeatureCollection(os.path.join(directory, 'Divergent.gpml'))
            self.assertGreater(len(divergent), 0)

    def test_a_declared_trench_exports_with_its_polarity(self):
        project = Project('test', *four_plates(), start_age_ma=1000., auto_trenches=False)
        history = History(project)
        history.commit({1: EAST}, events=[dict(kind='subduction', plates=[1, 2], overriding=2, points=STRETCH)])
        with tempfile.TemporaryDirectory() as directory:
            gplates.export(history, directory)
            trenches = [f for f in pygplates.FeatureCollection(os.path.join(directory, 'Convergent.gpml'))
                        if f.get_feature_type() == pygplates.FeatureType.gpml_subduction_zone]
            self.assertEqual(len(trenches), 1)
            trench = trenches[0]
            # Carried by the overriding plate, with the overriding plate on the side its polarity names.
            self.assertEqual(trench.get_reconstruction_plate_id(), 2)
            polarity = trench.get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
            line = trench.get_geometry()
            points = list(line.get_points())
            middle = len(points) // 2
            a, b = [np.array(p.to_xyz()) for p in points[middle:middle + 2]]   # away from the junctions
            left_of_line = sphere.unit(sphere.unit(a + b) + .05 * sphere.unit(np.cross(a, b)))
            lon, lat = sphere.to_lonlat(left_of_line)
            self.assertEqual(polarity, 'Right' if float(lon) < 0. else 'Left')   # plate 2 lies east of 0


if __name__ == '__main__':
    unittest.main()
