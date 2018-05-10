
# DIC using FFT-approach and optional pre- and post-processing

for 8bit, equal dim, single band images.
Pre-processing: Wallis filter.
Post-processing: RSME threshold, mean, median, vector filter.

A. Manconi & V. T. Bickel, 2.5.18.
andrea.manconi@erdw.ethz.ch / valentin.bickel@erdw.ethz.ch.
ETH Zurich / MPS Goettingen.

##### MIT License
##### Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi

Please cite this routine as:
#### Bickel, V.T., Manconi, A., Amann, F. (2018). Quantitative assessment of Digital Image Correlation methods to detect and monitor surface displacements of large slope instabilities. Remote Sensing Journal 2018.
________________________________________________________________________________________________________

### Required Matlab toolboxes:
- MATLAB R2017a
- Image Processing Toolbox
- (Communication System Toolbox for Mean and Median filter, optional)


### External/integrated functions/algorithms:
- dftregistration.m   by Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
Parts of their code has been taken from: J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 (1990).
Online: https://de.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation (Status: May 5th, 2018)
License included in "DIC/dft_registration_license.txt"

- wallis_filter.m,   underlying principle taken from Baltsavias, E. P. "Multiphoto geometrically constrained matching." Doctoral Thesis, ETH Zurich, (1991).


### Hardware requirements:

- Input image size of 30 MB requires at least ~5 GB RAM


## Instructions

1. Clone and unzip DIC_FFT_ETHZ_Repo

2. Place Master and Slave images in the "Input" folder

3. Customize setup for Preprocessing, DIC, and Postprocessing in run_pixel_offset.m:

          Inputs ----------------------------------------------------------------------------------
               - master = 'string' (name of master.type); # used for DIC
               - orig_m = 'string' (name of master.type); # used for plotting
               - slave = 'string' (name of slave.type);
               - orig_s = 'string' (name of slave.type);

          Preprocessing (Wallis filter & Co-registration) -----------------------------------------
               - wallis = 0/1  # wallis filter off/on
               - win = int. # window size, even numbers only
               - tarm = int. # target mean
               - tars = int. # target standard deviation
               - b = 0-1 # brightness enforcing constant
               - c = 0-1 # contrast enforcing constant

               - sp = int. # image split
               - co_os = int. # image oversampling

          DIC -------------------------------------------------------------------------------------
               - wi = 2^n, n = int. # window size in [pix]
               - os = int. # FFT oversampling factor
               - pix = float./int. # Ground sampling distance of used input imagery in [m]

          Postprocessing (RMSE, Mean, Vector, Median filter) --------------------------------------
               - filter = 1-4 # filter type, 1= RMSE threshold-, 2= Mean-, 3= Vector-, 4= Median filter

               - thr = 0-1 # threshold, 1 = no mask

               - mfws = [int. int.] # dimensions of mean filter window in [pix]
               - cut = int. # cut off value, OPTIONAL in [pix]
 
               - magcap = wi/int. or int. # tolerance-diff, values greater as value are cut, in [pix]
               - xcap = wi/int. or int.
               - ycap = wi/int. or int.

               - med = [int. int.] # dimensions of median filter window [pixel]

           Additional settings --------------------------------------------------------------------
               - scalax = [int. int.] # min and max values for the X displacement colorscale, in [m]
               - scalay = [int. int.] # min and max values for the Y displacement colorscale, in [m]

               - coppia = 'string' # name of output ascii file

               - skip_x = wi/int. # search window skip in X
               - skip_y = wi/int. # search window skip in Y

- Execute run_pixel_offset.m

- Collect displacement matrix ascii from "Output" folder

### Possible results

![](https://lh3.googleusercontent.com/WPJUSL6LNW25MJ8lrCsyz4GUkjhZgGNGDcMHltAd2iqR0rWJZSq8fh6IcwjxLaJos9U3QPqYH0-Q "DIC Example")
##### Application of various DIC algorithms to a landslide in the Swiss Alps; the magnitude of the displacement is displayed in meters. The central column represents the results derived with the available FFT-based approach. 
###### Taken from Bickel, V.T., Manconi, A., Amann, F. (2018). Quantitative assessment of Digital Image Correlation methods to detect and monitor surface displacements of large slope instabilities. Remote Sensing Journal 2018.

------------------
##### MIT License
##### Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
