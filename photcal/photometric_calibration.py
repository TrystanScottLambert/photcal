"""
Photometrically calibrating image by performing a reduction of the chi2 value for 
various color terms from a matching catalog and instrumental magnitudes from the host
telescope.
"""

from dataclasses import dataclass
from functools import cached_property
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats

from k_constant import calculate_k_constant_mag


@dataclass
class FilterMag:
    """Representation of a color observation"""
    name: str
    array: np.ndarray[float]
    array_err: np.ndarray[float]

@dataclass
class Color:
    """
    Represents the difference between two given colors. 
    The value will always be color_1 - color_2.
    """
    name: str
    filter_1: FilterMag
    filter_2: FilterMag

    @property
    def value(self):
        """Calculates the color"""
        return self.filter_1.array - self.filter_2.array

def get_unique_filters(list_filters: list[FilterMag]) -> list[FilterMag]:
    """Returns a unique list of given filters."""
    unique_filters = []
    for _filter in list_filters:
        if _filter not in unique_filters:
            unique_filters.append(_filter)
    return unique_filters


class Transformation:
    """All the inputs required to run the calibraton."""

    def __init__(
            self, observed_mag: FilterMag, catalog_mag: FilterMag, colors: list[Color]):
        """
        'observed_mag' is the FilterMag object of the observed magnitudes. Likewise, 
        'catalog_mag' is the FilterMag object of the catalog magnitudes. 'colors' is the 
        list of Color objects which represent the color terms. 
        """
        self.observed_mag = observed_mag
        self.catalog_mag = catalog_mag
        self.colors = colors

    @cached_property
    def sigma_squared(self) -> np.ndarray:
        """Determines the error term, sigma squared, from the errors of the given magnitudes."""
        filters = self._get_all_filters()
        sigma_array = np.zeros(len(self.observed_mag.array_err))
        for _filter in filters:
            sigma_array += _filter.array_err**2
        return sigma_array

    @cached_property
    def color_terms(self) -> np.ndarray:
        """Returns the array of color terms."""
        val = [np.ones(len(self.colors[0].value))]
        for color in self.colors:
            val.append(color.value)
        return val

    @property
    def delta_m(self) -> np.ndarray:
        """
        The difference between the observed magnitudes and the catalog magnitudes.
        """
        return self.observed_mag.array - self.catalog_mag.array

    def _get_all_filters(self):
        """
        Returns all the filters which are being used in the calculation, 
        including from the colors.
        """
        filters = []
        filters.append(self.observed_mag)
        filters.append(self.catalog_mag)
        for color in self.colors:
            filters.append(color.filter_1)
            filters.append(color.filter_2)
        return get_unique_filters(filters)

    def _create_matrix(self) -> np.array:
        """
        Creates the A matrix associated with the color terms. 
        """
        a_matrix = []
        for color_term in self.color_terms:
            row = []
            for other_color_term in self.color_terms:
                row.append(np.sum((color_term * other_color_term)/self.sigma_squared))
            a_matrix.append(row)
        return np.array(a_matrix)

    def _create_vector(self) -> np.array:
        """
        Creates the b vector which is associated with the 
        """
        b_vector = []
        for color_term in self.color_terms:
            b_vector.append(
                np.sum(((self.delta_m) * color_term) / (self.sigma_squared))
            )
        return np.array(b_vector)

    @cached_property
    def constants(self) -> np.ndarray:
        """Calculates the values of the constants."""
        a_matrix = self._create_matrix()
        b_vector = self._create_vector()
        return np.linalg.inv(a_matrix).dot(b_vector)

    @property
    def zero_point(self) -> float:
        """Determines the zero point for converting instrumental mags into catalog mags."""
        return -self.constants[0]

    def get_ap_corr_zero_point(self, aperture_radius: float, seeing: float) -> float:
        """
        Returns the zero point which has been aperture corrected. If the user wishes to use
        this zero point later on, then they will have to add a 'k_mag' correction to this zpt.
        This can be done using the calculate_k_constant_mag function in the k_constant module.

        'aperture_radius' is the radius of the apertures used for determining the
        magnitudes and the seeing of the image. Both 'aperture_radius' and 'seeing' must be 
        given in the same units.
        """
        return self.zero_point - calculate_k_constant_mag(aperture_radius, seeing)


    def diagnose(self):
        """Plot several diagnostic plots as well as information regarding the fit."""
        obs_mag_cal = self.observed_mag.array + self.zero_point
        cat_mag_cal = self.catalog_mag.array + self.constants[1]*(self.color_terms[1]) + self.constants[2]*(self.color_terms[2])
        mag_diff = obs_mag_cal - cat_mag_cal

        mean, median, std = sigma_clipped_stats(mag_diff)
        print('mean:', mean, 'median:', median, 'standard deviation:', std)

        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
        axs[0].plot(cat_mag_cal, self.observed_mag.array, 'o')
        axs[0].set_xlabel('Predicted catalog expected magnitudes')
        axs[0].set_ylabel('Instrumental magnitudes')

        axs[1].plot(obs_mag_cal, cat_mag_cal, 'o')
        axs[1].set_xlabel('Calibrated instrumental magnitudes')
        axs[1].set_ylabel('Predicted catalog expected magnitudes')

        axs[2].plot(self.catalog_mag.array, self.observed_mag.array-self.catalog_mag.array, 'o')
        axs[2].set_xlabel('Catalog Magnitudes')
        axs[2].set_ylabel('Instrumental Magnitudes - Catalog Magnitudes')

        plt.subplots_adjust(wspace=0.4)
        plt.show()
