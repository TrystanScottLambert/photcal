"""
Photometrically calibrating image by performing a reduction of the chi2 value for 
various color terms from a matching catalog and instrumental magnitudes from the host
telescope.
"""

from dataclasses import dataclass
from functools import cached_property, partial
import numpy as np


def photometric_transformation(
        m_inst: np.ndarray, constants: np.ndarray, color_term_list: list[np.ndarray]
        ) -> callable:
    """
    General form of the transformation from instrumental
    magnitude to photometrically corrected mag.

    Each color term should have the same number of elements as 
    the m_inst term. If we are using two color terms, for example, then 
    we should have an array of instrumental magnitudes x_mag as well as 
    some color terms  [r_mag - i_mag, i_mag - z_mag] for example. 
    """
    constant_terms = np.zeros(len(color_term_list[0]))
    for i, color_term in enumerate(color_term_list):
        constant_terms -= constants[i] * color_term
    return m_inst + constant_terms

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


@dataclass
class Settings:
    """All the inputs required to run the calibraton."""
    observed_mag: FilterMag
    catalog_mag: FilterMag
    colors: list[Color]

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

def create_matrix(settings: Settings) -> np.array:
    """
    Creates the A matrix associated with the color terms. 
    """
    a_matrix = []
    for color_term in settings.color_terms:
        row = []
        for other_color_term in settings.color_terms:
            row.append(np.sum((color_term * other_color_term)/settings.sigma_squared))
        a_matrix.append(row)
    return np.array(a_matrix)

def create_vector(settings: Settings) -> np.array:
    """
    Creates the b vector which is associated with the 
    """
    b_vector = []
    for color_term in settings.color_terms:
        b_vector.append(
            np.sum(((settings.delta_m) * color_term) / (settings.sigma_squared))
        )
    return np.array(b_vector)

def calculate_constants(a_matrix:np.matrix, b_vector: np.ndarray) -> np.ndarray:
    """Calculates the values of the constants."""
    return np.linalg.inv(a_matrix).dot(b_vector)

def get_photometric_transformation(settings: Settings) -> callable:
    """Determines the transformation for converting instrumental mags into catalog mags."""
    a_matrix = create_matrix(settings)
    b_vector = create_vector(settings)
    constants = calculate_constants(a_matrix, b_vector)
    transformation = partial(photometric_transformation, constants = constants)
    return transformation
