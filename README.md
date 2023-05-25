# photcal

This package can be used to photometrically calibrate data using stars identified in an astronomy image, and their corresponding apparent magnitudes in photometric surveys (i.e., GAIA or PANSTARRS).

In order to use the package the user needs to have a list of instrumental magnitudes with their errors as well as the corresponding magnitudes in a catalog with uncertainties. It is also important to have different color terms. So the user will also need additional magnitudes in other filters from the catalog.


## Installing the package 

The package is pip installable.

```pip install photcal```

## How it works

A full derivation of the method is available in the examples folder. 

## Using the package.
###Simple zero point determination.

The first step to use the package is to load in the different magnitude values for the different stars. To do this we provide the class ```FilterMag``` which takes in a name, a numpy array of magnitudes for the different stars, and an array of their associated errors. For example, if we wanted to calibrate GMOS r-band data with PANSTARRS we would identify the GMOS r-band stars using a software like sextractor or DAOstarfinder, as well as identify these stars in the PANSTARRS catalog. We could then create ```FilterMag``` objects:

```python
from photcat.photometric_calibration import FilterMag, Color, Settings, get_photometric_transformation
r_band_obs = FilterMag('r_gmos', r_gmos, r_gmos_err)
r_band_cat = FilterMag('r_pan', r_pan, r_pan_err)
g_band_cat = FilterMag('g_pan', g_pan, g_pan_err)
i_band_cat = FilterMag('i_pan', i_pan, i_pan_err)
```
We can then create out color terms which are used for the calibration. To do this we make use of the class ```Color``` which takes a name, a ```FilterMag``` object for the first magnitude in the catalog and a ```FilterMag``` object for the second magnitude in the catalog.:

```python
color_gr = Color('color_gr', g_band_cat, r_band_cat)
color_ir = Color('color_ir', i_band_cat, r_band_cat)
```

You can add as many color terms as you wish. We've used two color terms in this example.

All of these are then wrapped into one final class ```Transformation```  which takes the ```FilterMag``` of the image you want to calibrate (in this example GMOS r-band data), the ```FilterMag``` object of the magnitudes you want to calibrate to (in this case the PANSTARRS r-band data), and a list of ```Color``` objects which are all the color terms that you want to use.

```python
transformation = Transformation(r_band_obs, r_band_cat, [color_gr, color_ir])
```
At this point the calibration is finished and the user can use the zero point determination to calibrate their data.
```python
zpt = transformation.zero_point
r_gmos_calibrated  = r_gmos + zpt
```
It is also possible to check some diagnostic plots to make sure that the determined zpt value makes sense.
```python
transformation.diagnose()
```
### Taking into account aperture correction.
It is possible to calibrate the magnitudes whilst taking into account the size of the apertures which were used to determine the calibration. This may be useful if you need to change the aperture sizes for your source extraction later on and what to take this effect into account. 

To do this we have provided the ```k_constant```. All the steps of the calibration are the same except that now the zero point will have to include the radius of the apertures that were used as well as the seeing of your image. 

In our example we have a r-bad GMOS image with 0.8" seeing and used 2" apertures. These arguments must be passed with the same units. We can determine the zero point like

```python
zero_point_prime = transformation.get_ap_corr_zpt(aperture_radius=2, seeing=0.8)
```
Later, if we now decide to identify sources with 5" apertures instead, then we can calibrate the magnitudes like
```python
from photcal.k_constant import calculate_k_constant_mag
k_offset = calculate_k_constant_mag(aperture_radius = 5, seeing = 0.2)
zero_point = zero_point_prime + k_offset
r_gmos_calibrated  = r_gmos + zpt
```
If you won't change apertures then this process is unnecessary, however there are some selection criteria where this could be useful.
