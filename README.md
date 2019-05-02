
# Digital Image Correlation using an FFT-approach

#### State of the Art, validated, and calibrated DIC tool - for 8bit, equal dim, single- and multi-channel images, with on-demand geotiff information forwarding

#### Quantify displacement of features in a series of images over time & retrieve displacement velocities. Output data include CVP 2D-displacement maps, the displacement resultant (magnitude), displacement vectors, as well as a geo-rectified GIS-ready map (optional).
#### Built-in pre-processing routines:
(1) Wallis Filter: Dynamic Contrast Enhancement of both input images; improves the quality of the FFT correlation significantly and helps to suppress noise.

(2) Sub-Pixel Co-Registration: Aligns both input images on a sub-pixel level; significantly improves the quality of the correlation and the absolute accuracy of the displacement measurement.

#### Built-in post-processing routines:
(1) RMSE Threshold Filter: Filters DIC output based on a RMSE threshold.

(2) Mean Filter: Filters DIC output with a algorithmic mean kernel.

(3) Median Filter: Filters DIC output with a median kernel.

(4) Vector Filter: Filters DIC output based on the pixel neighborhood's displacement vector direction & magnitude.

<img src="https://github.com/bickelmps/DIC_FFT_ETHZ/blob/master/Figures/glacier.gif?raw=true">

-----------------------------------

V. T. Bickel & A. Manconi, May 2nd 2018

[valentin.bickel@erdw.ethz.ch / andrea.manconi@erdw.ethz.ch]

ETH Zurich / MPS Goettingen

##### MIT License - Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
Please cite this routine as:
#### Bickel, V.T.; Manconi, A.; Amann, F. Quantitative Assessment of Digital Image Correlation Methods to Detect and Monitor Surface Displacements of Large Slope Instabilities. Remote Sens. 2018, 10, 865.
http://www.mdpi.com/2072-4292/10/6/865
________________________________________________________________________________________________________

### Tool is being used by:
- Engineering Geology Group, ETH Zurich, Switzerland
- Swiss Seismological Service, ETH Zurich, Switzerland
- Department Planets and Comets, Max Planck Institute for Solar System Research, Germany
- [...]


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

2. Place Master and Slave image(s) in the "Input" folder (or copy images from 'Example' folder for a demo)

3. Customize setup for Preprocessing, DIC, and Postprocessing in run_pixel_offset.m, using the following variables. More information about the variables further below:

          Inputs ----------------------------------------------------------------------------------
               - geotiff = 0/1
               - epsg = int.
               - PCS = 'string'
               - master = 'string'
               - orig_m = 'string'
               - slave = 'string'
               - orig_s = 'string'
               - inputfilename = 'string'
               - outfilename = 'string'

          Preprocessing (Wallis filter & Co-registration) -----------------------------------------
               - wallis = 0/1
               - win = int.
               - tarm = int.
               - tars = int.
               - b = 0-1
               - c = 0-1

               - sp = int.
               - co_os = int.

          DIC -------------------------------------------------------------------------------------
               - wi = 2^n, n = int.
               - os = int.
               - pix = float./int.

          Postprocessing (RMSE, Mean, Vector, Median filter) --------------------------------------
               - filter = 1-4

               - thr = 0-1

               - mfws = [int. int.]
               - cut = int.
 
               - magcap = wi/int. or int.
               - xcap = wi/int. or int.
               - ycap = wi/int. or int.

               - med = [int. int.]

           Additional settings --------------------------------------------------------------------
               - scalax = [int. int.]
               - scalay = [int. int.]

               - coppia = 'string'

               - skip_x = wi/int.
               - skip_y = wi/int.

- Execute run_pixel_offset.m

- Collect displacement matrix ascii and deformation magnitude map tif from "Output" folder for further utilization of results (e.g. GIS)

## Input Parameter Description

          Inputs ----------------------------------------------------------------------------------
               - geotiff = indicate here, if you would like to forward geotiff information (=1) or not (=0)
               - epsg = in case you would like to forward geotiff information, enter the desired EPSG code here - you can get this number here: https://spatialreference.org/ref/epsg/
               - PCS = 
               - master = 
               - orig_m = 
               - slave = 
               - orig_s = 
               - inputfilename = 
               - outfilename = 

          Preprocessing (Wallis filter & Co-registration) -----------------------------------------
               - wallis = 
               - win = 
               - tarm = 
               - tars = 
               - b = 
               - c = 

               - sp = 
               - co_os = 

          DIC -------------------------------------------------------------------------------------
               - wi = 
               - os = 
               - pix = 

          Postprocessing (RMSE, Mean, Vector, Median filter) --------------------------------------
               - filter = 

               - thr = 

               - mfws = 
               - cut = 
 
               - magcap = 
               - xcap = 
               - ycap = 

               - med = 

           Additional settings --------------------------------------------------------------------
               - scalax = 
               - scalay = 

               - coppia = 

               - skip_x = 
               - skip_y = 

### Possible results

##### Advancing glacier, one image per day:

<img src="https://github.com/bickelmps/DIC_FFT_ETHZ/blob/master/Figures/glacier.gif?raw=true">

------------------
##### MIT License - Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
