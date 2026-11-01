"""The GPlates-to-goSPL pipeline on a small project made here: one continent on one plate."""
import json
import os
import tempfile
import unittest

import numpy as np
import pygplates

from gplates_gospl import build, motion
from gplates_gospl.elevation import Mesh, elevation, ocean_depth_m
from gplates_gospl.mesh import icosphere
from gplates_gospl.project import Project

ROLES = {'ocean': ['ocean.gpml'], 'continent': ['continents.gpml'], 'orogeny_active': ['orogenies.gpml']}


def make_project(directory):
    """Plate 101 turns east 1 degree per Myr about the north pole, from 100 Ma (where it is drawn) to
    the present; it carries a continent with a range down its middle, beside an ocean of crust
    formed at 150 Ma on the fixed plate 102."""
    with open(os.path.join(directory, 'rotation.rot'), 'w') as f:
        f.write('101 0.0 90.0 0.0 100.0 000 !\n101 100.0 90.0 0.0 0.0 000 !\n')
        f.write('102 0.0 90.0 0.0 0.0 000 !\n102 100.0 90.0 0.0 0.0 000 !\n')

    def save(name, kind, plate, ring, begin=1000.):
        feature = pygplates.Feature(pygplates.FeatureType.create_gpml(kind))
        feature.set_geometry(pygplates.PolygonOnSphere(ring))
        feature.set_reconstruction_plate_id(plate)
        feature.set_valid_time(begin, -1e-3)
        pygplates.FeatureCollection([feature]).write(os.path.join(directory, name))

    save('continents.gpml', 'ContinentalCrust', 101, [(-20, -20), (-20, 20), (20, 20), (20, -20)])
    save('orogenies.gpml', 'OrogenicBelt', 101, [(-15, -3), (-15, 3), (15, 3), (15, -3)])
    save('ocean.gpml', 'OceanicCrust', 102, [(-60, 120), (60, 120), (60, -120), (-60, -120)], begin=150.)


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.mkdtemp()
        make_project(cls.directory)
        cls.project = Project(cls.directory, ROLES)
        nodes, cls.tris = icosphere(5)
        cls.mesh = Mesh(nodes)

    def test_the_icosphere_is_closed(self):
        nodes, tris = icosphere(3)
        self.assertEqual(len(nodes) - 3 * len(tris) // 2 + len(tris), 2)   # Euler characteristic
        normals = np.cross(nodes[tris[:, 1]] - nodes[tris[:, 0]], nodes[tris[:, 2]] - nodes[tris[:, 0]])
        self.assertTrue(np.all(np.einsum('ij,ij->i', normals, nodes[tris[:, 0]]) > 0.))

    def test_land_stands_up_ranges_rise_and_the_ocean_deepens_with_age(self):
        z = elevation(self.project, 100., self.mesh)['z']
        lon, lat = self.mesh.lon, self.mesh.lat
        interior = (np.abs(lat) < 12.) & (np.abs(lon - 10.) < 4.)
        self.assertTrue(np.all(z[interior] > 0.))
        ranges = (np.abs(lat) < 12.) & (np.abs(lon) < 2.)
        self.assertGreater(z[ranges].mean(), z[interior].mean() + 1000.)
        far = (np.abs(lat) < 40.) & (np.abs(((lon + 360.) % 360.) - 180.) < 30.)
        np.testing.assert_allclose(z[far], ocean_depth_m(50.), atol=1.)   # crust of 150 Ma, at 100 Ma

    def test_nodes_go_with_the_plate_of_their_crust(self):
        ids = motion.plate_ids(self.project, 100., self.mesh)
        on_land = (np.abs(self.mesh.lat) < 15.) & (np.abs(self.mesh.lon) < 15.)
        self.assertTrue(np.all(ids[on_land] == 101))
        moved = motion.carried(self.project, ids, self.mesh.xyz, 100., 90.)
        dlon = (np.degrees(np.arctan2(moved[:, 1], moved[:, 0])) - self.mesh.lon + 540.) % 360. - 180.
        np.testing.assert_allclose(dlon[on_land & (np.abs(self.mesh.lat) < 80.)], 10., atol=1e-6)
        self.assertTrue(np.allclose(moved[ids == 102], self.mesh.xyz[ids == 102]))

    def test_the_land_drains_to_the_sea(self):
        from gplates_gospl.elevation import drained
        nodes, tris = icosphere(5)
        mesh = Mesh(nodes, tris)
        lat = np.degrees(np.arcsin(nodes[:, 2]))
        z = np.where(lat > 0., 500. - 400. * np.exp(-((lat - 45.) / 10.) ** 2), -3000.)   # a basin inland
        out = drained(z, mesh.neighbours)
        self.assertTrue(np.all(out >= z))
        self.assertTrue(np.all(out[z <= 0.] == z[z <= 0.]))
        # Every land node has a neighbour lower than itself, or the sea: no closed depression is left.
        n = mesh.neighbours
        for node in np.where(out > 0.)[0]:
            near = n.indices[n.indptr[node]:n.indptr[node + 1]]
            self.assertTrue((out[near] < out[node]).any(), node)

    def test_smoothing_keeps_land_and_sea_apart(self):
        from gplates_gospl.fields import smoothed
        nodes, _ = icosphere(4)
        land = nodes[:, 0] > 0.
        values = np.where(land, 100., -100.) + np.sin(60. * nodes[:, 2])
        out = smoothed(values, nodes, 800., within=land)
        self.assertTrue(np.all(out[land] > 98.) and np.all(out[~land] < -98.))   # no blur across the coast
        self.assertLess(np.std(out[land]), np.std(values[land]) / 3.)              # the ripple is smoothed

    def test_the_forcing_carries_goSPLs_own_surface_to_each_map(self):
        with tempfile.TemporaryDirectory() as out:
            ages, maps = build.build(self.directory, out, older=100., younger=70., every=10., level=5, roles=ROLES)
            data = np.load(os.path.join(out, 'input', 'maps.npz'))
            mesh = np.load(os.path.join(out, 'input', 'mesh.npz'))
            nodes = mesh['v'] / build.RADIUS_M
            # goSPL's fast mode, step by step: lift over the interval, then carry by interpolation.
            surface = mesh['z']
            for k in range(len(ages) - 1):
                forcing = np.load(os.path.join(out, 'input', 'forcing_%03d.npz' % k))
                years = (ages[k] - ages[k + 1]) * 1e6
                moved = mesh['v'] + forcing['hdisp'] * years
                surface = build.advected(surface + forcing['t'] * years, mesh['v'], moved)
                np.testing.assert_allclose(surface, data['expected'][k + 1], atol=1e-6)
            # Carried onto the fixed mesh the surface softens by about a node, no more.
            self.assertLess(np.mean((surface > 0) != (data['z'][-1] > 0)), .05)
            with open(os.path.join(out, 'forcing-check.yml')) as f:
                config = json.load(f)
            self.assertTrue(config['domain']['fast'])
            self.assertEqual(config['time']['end'], 30e6)
            self.assertEqual([e['start'] for e in config['tectonics']], [0., 10e6, 20e6])
            self.assertTrue(np.allclose(np.linalg.norm(nodes, axis=1), 1.))


class Rocke3dTests(unittest.TestCase):
    def test_a_climate_goes_from_the_extractors_csv_onto_the_mesh(self):
        from gplates_gospl import rocke3d
        nodes, _ = icosphere(5)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'climate.csv')
            rocke3d.random_climate_csv(path, seed=3)
            with open(path) as f:
                header = f.readline().strip().split(',')
            self.assertEqual(header[:4], ['j', 'i', 't01_c', 'p01_mm_day'])
            self.assertEqual(len(header), 29)
            grid = rocke3d.read_precipitation(path)
        self.assertEqual(grid.shape, (46, 72))
        self.assertTrue(np.all(grid[0] == grid[0, 0]) and np.all(grid[-1] == grid[-1, 0]))   # caps
        rain = rocke3d.rain_on_mesh(grid, nodes)
        self.assertTrue(np.all(rain > 0.))
        self.assertTrue(.3 < rain.mean() < 3.)   # m/yr
        lat = np.degrees(np.arcsin(nodes[:, 2]))
        self.assertGreater(rain[np.abs(lat) < 8.].mean(), rain[np.abs(np.abs(lat) - 25.) < 5.].mean())

    def test_the_month_lengths_weight_the_annual_mean(self):
        from gplates_gospl import rocke3d
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'climate.csv')
            with open(path, 'w') as f:
                f.write(','.join(['j', 'i'] + ['t%02d_c,p%02d_mm_day' % (m, m) for m in range(1, 13)]
                                 + ['u_face_cm_s', 'v_face_cm_s', 'w_top_cm_s']) + '\n')
                for j in range(46):
                    for i in range(72):
                        f.write('%d,%d,' % (j, i) + ','.join('0,%g' % (10. if m == 0 else 0.) for m in range(12))
                                + ',nan,nan,nan\n')
            grid = rocke3d.read_precipitation(path, month_days=[30] + [10] * 11)
        np.testing.assert_allclose(grid, 10. * 30. / 140. * 365.25 / 1000.)

    def test_a_world_goes_onto_aeries_grid(self):
        from gplates_gospl import rocke3d
        nodes, _ = icosphere(6)
        z = np.where(nodes[:, 2] > .5, 1000., -4000.)   # land north of 30 N
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'grid.npz')
            rocke3d.grid_npz(path, nodes, z)
            grid = np.load(path)
            self.assertEqual(grid['land_fraction'].shape, (46, 72))
            self.assertAlmostEqual(float(grid['cell_area_steradian'][1:-1].sum() + grid['cell_area_steradian'][0, 0]
                                         + grid['cell_area_steradian'][-1, 0]), 4. * np.pi, places=6)
            lat = grid['latitude_deg_n']
            self.assertTrue(np.all(grid['land_fraction'][lat > 34.] == 1.))
            self.assertTrue(np.all(grid['land_fraction'][lat < 26.] == 0.))
            np.testing.assert_allclose(grid['dominant_surface_elevation_m'][lat > 34.], 1000.)


class RainShadowTests(unittest.TestCase):
    def test_a_range_wets_its_windward_side_and_dries_its_lee(self):
        from gplates_gospl import rain
        from gplates_gospl.fields import smoothed
        nodes, _ = icosphere(6)
        lon = np.degrees(np.arctan2(nodes[:, 1], nodes[:, 0]))
        lat = np.degrees(np.arcsin(nodes[:, 2]))
        # A north-south range at 0 E, in the westerlies (30-60 N), on land from 20 W to 20 E.
        z = np.where((np.abs(lon) < 20.) & (lat > 20.) & (lat < 70.), 200. + 3000. * np.exp(-(lon / 3.) ** 2), -4000.)
        base = np.ones(len(nodes))
        shaped = rain.orographic(base, z, nodes)
        band = (lat > 40.) & (lat < 50.)
        windward = band & (lon > -4.) & (lon < 1.5)     # the westerlies come from the west
        lee = band & (lon > 1.5) & (lon < 8.)
        self.assertGreater(shaped[windward].mean(), 1.1)
        self.assertLess(shaped[lee].mean(), .8)
        np.testing.assert_allclose(smoothed(shaped, nodes, 1500.)[band], 1., atol=.15)   # moved, not made

    def test_the_trades_blow_from_the_east_and_the_westerlies_from_the_west(self):
        from gplates_gospl import rain
        nodes = np.array([[1., 0., 0.], [np.cos(np.radians(45.)), 0., np.sin(np.radians(45.))],
                          [np.cos(np.radians(30.)), 0., np.sin(np.radians(30.))]])
        wind = rain.winds(nodes)
        self.assertLess(wind[0, 1], -.99)      # at 0 N, 0 E, blowing west
        self.assertGreater(wind[1, 1], .99)    # at 45 N, blowing east
        self.assertLess(np.linalg.norm(wind[2]), 1e-9)   # calm at 30 N


if __name__ == '__main__':
    unittest.main()
