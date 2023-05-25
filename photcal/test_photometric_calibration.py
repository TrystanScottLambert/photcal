"""
Tests for photometric calibration
"""

import sys
import unittest
import numpy as np
from astropy.table import Table

from photometric_calibration import FilterMag, Transformation, Color, get_unique_filters

def read_in_test_data():
    """Reads in the test file."""
    i_obs, i_obs_err, i_cat, i_cat_err, z_cat, z_cat_err, r_cat, r_cat_err =\
            np.loadtxt(
        'test_cat.dat', unpack=True, usecols=(2, 3, 12, 13, 14, 15, 10, 11)
    )
    return i_obs, i_obs_err, i_cat, i_cat_err, z_cat, z_cat_err, r_cat, r_cat_err


class TestFilterMag(unittest.TestCase):
    """Testing that the class 'FilterMag' is reading in correctly."""

    def test_reading_in(self):
        """testing that filter mag is reading in correctly."""
        data = read_in_test_data()
        filter_mag = FilterMag('test',data[0], data[1])
        self.assertIsInstance(filter_mag, FilterMag)
        self.assertEqual(filter_mag.name, 'test')
        self.assertListEqual(list(filter_mag.array), list(data[0]))
        self.assertListEqual(list(filter_mag.array_err), list(data[1]))

class TestColor(unittest.TestCase):
    """Testing that the class 'Color' is reading in correctly."""
    def test_reading_in(self):
        """testing that the test color is reading in correctly and determining correct values."""
        data = read_in_test_data()
        i_band_cat = FilterMag('i_cat', data[1], data[2])
        z_band_cat = FilterMag('z_cat', data[4], data[5])
        zi_color = Color('zi_color', i_band_cat, z_band_cat)
        self.assertListEqual(list(zi_color.value), list(i_band_cat.array - z_band_cat.array))

class TestTransformation(unittest.TestCase):
    """Testing that the Transformation class is working correctly."""

    def create_transformation_object(self):
        """Creates a Transformation object which can be used in tests."""
        data = read_in_test_data()
        i_band_obs = FilterMag('i_obs', data[0], data[1])
        i_band_cat = FilterMag('i_cat', data[2], data[3])
        z_band_cat = FilterMag('z_cat', data[4], data[5])
        r_band_cat = FilterMag('r_cat', data[6], data[7])
        iz_color = Color('iz_color', i_band_cat, z_band_cat)
        ri_color = Color('ir_color', r_band_cat, z_band_cat)
        transformation = Transformation(i_band_obs, i_band_cat, [iz_color, ri_color])
        return transformation

    def setup_test(self):
        """Creates the data and the Transformation object which are needed for tests."""
        data = read_in_test_data()
        transformation = self.create_transformation_object()
        return data, transformation

    def test_sigma_squared(self):
        """testing if the error term is being calculated correctly"""
        data, transformation = self.setup_test()
        self.assertListEqual(list(transformation.sigma_squared),
                             list(data[1]**2 + data[3]**2 + data[5]**2 + data[7]**2))

    def test_delta_m(self):
        """Testing that delta m give the correct value."""
        data, transformation = self.setup_test()
        self.assertListEqual(list(transformation.delta_m), list(data[0] - data[2]))

    def test_color_terms(self):
        """testing that the color terms are correct"""
        data = read_in_test_data()
        i_band_obs = FilterMag('i_obs', data[0], data[1])
        i_band_cat = FilterMag('i_cat', data[2], data[3])
        z_band_cat = FilterMag('z_cat', data[4], data[5])
        r_band_cat = FilterMag('r_cat', data[6], data[7])
        iz_color = Color('iz_color', i_band_cat, z_band_cat)
        ri_color = Color('ir_color', r_band_cat, z_band_cat)
        transformation = Transformation(i_band_obs, i_band_cat, [iz_color, ri_color])
        color_terms = [
            np.ones(len(i_band_cat.array)),
            i_band_cat.array - z_band_cat.array,
            r_band_cat.array - z_band_cat.array]
        for i, color_term in enumerate(color_terms):
            self.assertListEqual(list(transformation.color_terms[i]), list(color_term))
    
    def test_matrix(self):
        """testing matrix works."""
        _, transformation = self.setup_test()

        a_matrix = transformation._create_matrix()

        self.assertEqual(a_matrix.shape, (3, 3))
        self.assertEqual(a_matrix[1,1], np.sum((transformation.color_terms[1]*transformation.color_terms[1])/transformation.sigma_squared))
        self.assertEqual(a_matrix[1,2], np.sum((transformation.color_terms[1]*transformation.color_terms[2])/transformation.sigma_squared))
        self.assertEqual(a_matrix[0,0], np.sum((transformation.color_terms[0]*transformation.color_terms[0])/transformation.sigma_squared))
        self.assertEqual(a_matrix[-1,1], np.sum((transformation.color_terms[-1]*transformation.color_terms[1])/transformation.sigma_squared))
        self.assertEqual(a_matrix[-1,-1], np.sum((transformation.color_terms[-1]*transformation.color_terms[-1])/transformation.sigma_squared))
    
    def test_b_vector(self):
        """testing that the b vector is correct."""
        _, transformation = self.setup_test()
        b_vector = transformation._create_vector()

        for i, color_term in enumerate(transformation.color_terms):
            self.assertEqual(b_vector[i], np.sum(((transformation.delta_m) * (color_term))/transformation.sigma_squared))


class TestUniqueFilters(unittest.TestCase):
    "Testing that unique_filters returns a unique list of given FilterMags"

    def test_unique(self):
        """Doing the test."""
        data = read_in_test_data()
        i_band_obs = FilterMag('i_obs', data[0], data[1])
        i_band_cat = FilterMag('i_cat', data[1], data[2])
        z_band_cat = FilterMag('z_cat', data[4], data[5])
        r_band_cat = FilterMag('r_cat', data[6], data[7])
        filter_list = [i_band_cat, i_band_obs, i_band_cat, z_band_cat, z_band_cat, r_band_cat]
        unique_list = get_unique_filters(filter_list)
        self.assertListEqual(unique_list, [i_band_cat, i_band_obs, z_band_cat, r_band_cat])


class TestAgainstExample(unittest.TestCase):
    """We use an example which was determined by hand to test that the code is behaving in a similar way"""

    #manual method
    gmos_panss = Table.read('test_cat.dat', format='ascii')

    r_gmos = gmos_panss['rmag']
    r_pan = gmos_panss['rApMag']
    g_pan = gmos_panss['gApMag']
    i_pan = gmos_panss['iApMag']

    re_gmos = gmos_panss['rmagerr']
    re_pan = gmos_panss['rApMagErr']
    ge_pan = gmos_panss['gApMagErr']
    ie_pan = gmos_panss['iApMagErr']
    sigma = re_gmos**2 + re_pan**2 + ge_pan**2 + ie_pan**2

    A11 = np.sum(1/sigma)
    A12 = np.sum((g_pan-r_pan)/sigma)
    A13 = np.sum((i_pan-r_pan)/sigma)
    A21 = np.sum((g_pan-r_pan)/sigma)
    A22 = np.sum((g_pan-r_pan)**2/sigma)
    A23 = np.sum((i_pan-r_pan)*((g_pan-r_pan)/sigma))
    A31 = np.sum((i_pan-r_pan)/sigma)
    A32 = np.sum((i_pan-r_pan)*((g_pan-r_pan)/sigma))
    A33 = np.sum((i_pan-r_pan)**2/sigma)

    A = np.array([[A11, A12, A13],
                [A21, A22, A23],
                [A31, A32, A33]])

    b1 = np.sum((r_gmos-r_pan)/sigma)
    b2 = np.sum((r_gmos-r_pan)*((g_pan-r_pan)/sigma))
    b3 = np.sum((r_gmos-r_pan)*((i_pan-r_pan)/sigma))

    B = np.array([b1,b2,b3])

    const = np.linalg.inv(A).dot(B)
    C1 = float(const[0])
    C2 = float(const[1])
    C3 = float(const[2])

    cal_gmos_r =r_gmos - C1

    #New method
    r_band_obs = FilterMag('r_gmos', np.array(r_gmos), np.array(re_gmos))
    r_band_cat = FilterMag('r_pan', np.array(r_pan), np.array(re_pan))
    g_band_cat = FilterMag('g_pan', np.array(g_pan), np.array(ge_pan))
    i_band_cat = FilterMag('i_pan', np.array(i_pan), np.array(ie_pan))
    color_gr = Color('color_gr', g_band_cat, r_band_cat)
    color_ir = Color('color_ir', i_band_cat, r_band_cat)
    transformation = Transformation(r_band_obs, r_band_cat, [color_gr, color_ir])


    def test_sigma_squared(self):
        """testing that both values of sigma squared are the same."""
        self.assertListEqual(list(self.transformation.sigma_squared), list(self.sigma))

    def test_matrix(self):
        """testing that the matrix values are the same from the manual and generalized example."""

        a_matrix = self.transformation._create_matrix()
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(self.A[i, j], a_matrix[i,j])

    def test_vector(self):
        """Testing that the vector is the same as the manual method"""

        b_vector = self.transformation._create_vector()
        for i, value in enumerate(self.B):
            self.assertAlmostEqual(value, b_vector[i])

    def test_constants(self):
        """Testing that the constants are the same as the given example."""
        constant_values = self.transformation.constants
        for i, value in enumerate(self.const):
            self.assertAlmostEqual(value, constant_values[i])

    def test_transform(self):
        """Return the transform and determine if it is correct."""
        corrected_mags = np.array(self.r_gmos) + self.transformation.zero_point
        for i, value in enumerate(self.cal_gmos_r):
            self.assertAlmostEqual(value, corrected_mags[i])
        

if __name__ == '__main__':
    unittest.main()
