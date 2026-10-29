"""Elevation from the world's features, and the maps made from it."""
import math
import os
import tempfile
import unittest

import numpy as np
from scipy.io import netcdf_file

from paradise import relief
from paradise.core import sphere
from paradise.history import History, Project
from tests.test_paradise import four_plates, grid_lonlat


def with_a_continent():
    plates, crust = four_plates()
    lon, lat = grid_lonlat()
    crust[sphere.arc_km(sphere.from_lonlat(lon, lat), sphere.from_lonlat(-90., 0.)[None, None]) < 2500.] = 1
    return Project('test', plates, crust, start_age_ma=1000., ocean_age_myr=60.)


class ElevationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.history = History(with_a_continent())
        cls.z, cls.lon, cls.lat = relief.elevation(cls.history.current, width=360)

    def at(self, lon, lat):
        row = int((90. - lat) / 180. * self.z.shape[0])
        col = int((lon + 180.) / 360. * self.z.shape[1])
        return self.z[row, col]

    def test_continent_stands_above_the_sea_and_old_ocean_lies_deep(self):
        self.assertGreater(self.at(-90., 0.), 0.)
        expected = -1000. * (relief.RIDGE_KM + relief.DEPTH_PER_ROOT_MYR * math.sqrt(60.))
        self.assertAlmostEqual(self.at(90., 0.), expected, delta=500.)

    def test_land_share_is_about_the_continent(self):
        weight = np.cos(np.radians(self.lat))
        land = float((weight * (self.z > 0)).sum() / weight.sum())
        cap = (1. - math.cos(2500. / sphere.RADIUS_KM)) / 2.
        self.assertAlmostEqual(land, cap, delta=.3 * cap)


class AerieTests(unittest.TestCase):
    def test_the_grid_follows_aeries_contract_and_the_topo_file_its_identity(self):
        history = History(with_a_continent())
        with tempfile.TemporaryDirectory() as directory:
            files = relief.write_maps(history, directory, width=360)
            names = {os.path.basename(f) for f in files}
            self.assertEqual(names, {'elevation.png', 'heightmap-16bit.png', 'rocke3d-4x5.npz', 'rocke3d-TOPO.nc',
                                     'maps.json'})
            grid = np.load(os.path.join(directory, 'rocke3d-4x5.npz'))
            for key in ('cell_area_steradian', 'land_fraction', 'mean_signed_elevation_m', 'mean_land_elevation_m',
                        'mean_ocean_elevation_m', 'dominant_surface_elevation_m'):
                self.assertEqual(grid[key].shape, (46, 72), key)
            self.assertEqual(grid['latitude_bounds_deg'].shape, (47,))
            self.assertAlmostEqual(grid['cell_area_steradian'].sum(), 4. * math.pi, places=9)
            self.assertTrue(np.all(grid['land_fraction'][0] == grid['land_fraction'][0, 0]))   # the cap, repeated
            self.assertTrue(np.all(np.isnan(grid['mean_land_elevation_m'][grid['land_fraction'] == 0.])))
            with netcdf_file(os.path.join(directory, 'rocke3d-TOPO.nc'), 'r', mmap=False) as f:
                v = {k: np.array(f.variables[k][:], float) for k in f.variables}
            np.testing.assert_allclose(v['focean'] + v['fgrnd'], 1., atol=1e-6)
            np.testing.assert_allclose(v['zsolid'], v['zatmo'] - v['focean'] * v['hocean'], atol=1e-2)
            self.assertGreaterEqual(v['hocean'][v['focean'] > 0].min(), 10.)


if __name__ == '__main__':
    unittest.main()
