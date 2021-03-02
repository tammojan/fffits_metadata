## Add FITS metadata to FF files

This is a small script to add metadata to fits files ("FF files") created by [RMS](https://github.com/CroatianMeteorNetwork/RMS). It searches for the most applicable platepar (calibration solutions) and approximates them by a FITS WCS, so that the FITS files can be read in any FITS reader while keeping coordinates.

Additional information added to the metadata:

 * Exact time (which is also in the filename).
 * Station location (latitude and longitude only, rounded to two decimals)
 * Reference to Global Meteor Network
 * Information about number of frames, framerate and exposure

Example metadata:
```
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                    8 / array data type
NAXIS   =                    0 / number of array dimensions
EXTEND  =                    T
NROWS   =                  720
NCOLS   =                 1280
NBITS   =                    8
NFRAMES =                  256
FIRST   =               234240
CAMNO   = 'NL000D  '
FPS     =                 25.0
OBSERVER= 'NL000D  '
INSTRUME= 'Global Meteor Network'
MJD-OBS =    59211.79409037037
DATE-OBS= '2020-12-28T19:03:29.408'
EXPTIME =                10.24
SITELONG=                 6.36
SITELAT =                52.84
```

### Usage
```

```usage: add_fits_metadata.py [-h] dir_path
add_fits_metadata.py: error: the following arguments are required: dir_path
```