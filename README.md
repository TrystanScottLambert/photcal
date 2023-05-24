# photcal

This package can be used to photometrically calibrate data using stars identified in an astronomy image, and their corresponding apparent magnitudes in photometric surveys (i.e., GAIA or PANSTARRS).

In order to use the package the user needs to have a list of instrumental magnitudes with their errors as well as the corresponding magnitudes in a catalog with uncertainties. It is also important to have different color terms. So the user will also need additional magnitudes in other filters from the catalog.


## Installing the package 

The package is pip installable.

```pip install photcal```

## How it works

A full derivation of the method is available in the examples folder. 

## Using the package. 

The first step to use the pacakge is to load in the different magnitude values for the different stars. To do this we provide the class ```FilterMag``` which takes in a name, a numpy array of magnitudes for the different stars, and an array of their associated errors. For example, if we wanted to calibrate GMOS r-band data with PANSTARRS we would identify the GMOS r-band stars using a software like sextractor or DAOstarfinder, as well as identify the stars in the PANSTARRS catalog. We could then create the ```FilterMag``` objects like:

```python
from photcat.photometric_calibration import FilterMag, Color, Settings, get_photometric_transformation
r_band_obs = FilterMag('r_gmos', r_gmos, np.array(re_gmos))
r_band_cat = FilterMag('r_pan', r_pan, np.array(re_pan))
g_band_cat = FilterMag('g_pan', g_pan, np.array(ge_pan))
i_band_cat = FilterMag('i_pan', i_pan, np.array(ie_pan))
```
We can then create out color terms which are used for the calibration. To do this we make use of the class ```Color``` which takes a name, a ```FilterMag``` object for the first magnitude and a ```FilterMag``` object for the second magnitude:

```python
color_gr = Color('color_gr', g_band_cat, r_band_cat)
color_ir = Color('color_ir', i_band_cat, r_band_cat)
```

You can add as many color terms as you wish. 

All these are then wrapped into one final class ```Settings```  which takes the ```FilterMag``` of the image you want to calibrate (in this example GMOS r-band data), the ```FilterMag``` object of the magnitudes you want to calibrate to (in this case the PANSTARRS r-band data), and a list of ```Color``` objects which are all the color terms that you want to use.

```python
settings = Settings(r_band_obs, r_band_cat, [color_gr, color_ir])
```

The transformation function from instrumental magnitudes to calibrated magnitudes can be determined using the ```get_photometric_transformation``` function.
```python
transform_inst_to_cal = get_photometric_transformation(settings)
```
```get_photometric_transformation``` returns a function which you can then use to conver either a single magnitude value, or indeed, an entire array:

```python
calibrated_r_magnitudes = transform_inst_to_cal(r_gmos)
```
