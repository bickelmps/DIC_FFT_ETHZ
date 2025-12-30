
# Digital Image Correlation using an FFT-approach

#### State of the Art, validated, and calibrated DIC tool - for 8bit, equal dim, single- and multi-channel images, with on-demand geotiff information forwarding

#### Quantify displacement of features (landslides, glaciers, etc.) in a series of images over time & retrieve displacement velocities and directions. Output data include CVP 2D-displacement maps, displacement vectors, as well as geo-rectified GIS-ready (.tif and .txt) output (optional).

V. T. Bickel & A. Manconi, May 22nd 2020 (v3) | ETH Zurich / MPS Goettingen

Please direct bugs/questions/feedback to: [valentin.bickel@erdw.ethz.ch / andrea.manconi@erdw.ethz.ch]
________________________________________________________________________________________________________
Please cite this routine as:
#### Bickel, V.T.; Manconi, A.; Amann, F. Quantitative Assessment of Digital Image Correlation Methods to Detect and Monitor Surface Displacements of Large Slope Instabilities. Remote Sens. 2018, 10, 865.
http://www.mdpi.com/2072-4292/10/6/865
##### MIT License - Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi

<img src="https://github.com/bickelmps/DIC_FFT_ETHZ/blob/master/Figures/example3.png?raw=true">

##### Figure. DIC_FFT application example for a slope instability in the Swiss Alps: (a) slope instability displacement field over 19 years returned using DIC_FFT and two optical images, plotted on a hillshade. Images taken by an ebee drone (fixed wing) and from a plane. Blank areas indicate portions of the slope with displacements below DIC accuracy or without information. All runs were performed with a scan window size of 64 × 64 pixels and a spatial vector filter. (b) Interpolation of the Cuolm da Vi slope displacement field based on 23 GNSS measurements as direct validation. Results shows good agreement with the results produced using DIC_FFT, while DIC_FFT produces a much higher spatial resolution: The fast-moving section visible in the West of the area failed in 2016; using GNSS data alone, the pre-failure displacement is not visible. Modified from Bickel et al., 2018. More examples at the bottom of the readme!

## Graphical Abstract
<img src="https://github.com/bickelmps/DIC_FFT_ETHZ/blob/master/Figures/graphic_abs.png?raw=true">
________________________________________________________________________________________________________

#### Built-in pre-processing routines:
(1) Wallis Filter: Dynamic Contrast Enhancement of both input images; improves the quality of the FFT correlation significantly and helps to suppress noise.

(2) Sub-Pixel Co-Registration: Aligns both input images on a sub-pixel level; significantly improves the quality of the correlation and the absolute accuracy of the displacement measurement.

#### Built-in post-processing routines:
(1) RMSE Threshold Filter: Filters DIC output based on a RMSE threshold.

(2) Mean Filter: Filters DIC output with a algorithmic mean kernel.

(3) Median Filter: Filters DIC output with a median kernel.

(4) Vector Filter: Filters DIC output based on the pixel neighborhood's displacement vector direction & magnitude.

<img src="https://github.com/bickelmps/DIC_FFT_ETHZ/blob/master/Figures/glacier.gif?raw=true">

##### Figure. Example of DIC_FFT application for glacier displacement monitoring. The subplots show 1) primary and 2) secondary images, 3) the displacement magnitude, as well as 4) the displacement vectors. More examples at the bottom of the readme!
________________________________________________________________________________________________________

### Tool is being used by:
- Engineering Geology Group, ETH Zurich, Switzerland
- Swiss Seismological Service, ETH Zurich, Switzerland
- WSL Institute for Snow and Avalanche Research SLF, Switzerland
- Department Planets and Comets, Max Planck Institute for Solar System Research, Germany
- Landslide Research Group, Technical University of Munich, Germany
- Department of Engineering Geology and Hydrogeology, RWTH Aachen, Germany
- School of Earth and Environment, University of Canterbury, New Zealand
- Centre for Hydrogeology and Geothermics, University of Neuchatel, Switzerland
- Institute of Photogrammetry and Geoinformation, Leibniz University Hannover, Germany
- German Research Centre for Geosciences, Helmholtz Centre Potsdam, Germany
- CSD Engineers, Switzerland
- Geology & Geophysics, University of Utah, United States of America
- State Key Laboratory for Tunnel Engineering, China University of Mining and Technology-Beijing, Beijing, China
- [...] 

### Publications

- Zheng, X., Tao, Z., Yu, H. Model Study on the Failure Mechanism of Low-Inclination Anti-dip Slopes with Weak Interlayers. Rock Mechanics and Rock Engineering, 2025 (https://doi.org/10.1007/s00603-025-05169-1)
- Louis, C., Halloran, L., Roques, C. Seasonal and diurnal freeze–thaw dynamics of a rock glacier and their impacts on mixing and solute transport. Hydrology and Earth System Sciences, 2025 (https://doi.org/10.5194/hess-29-1505-2025)
- Manconi, A., Jones, N., Loew, S., et al. Monitoring surface deformation with spaceborne radar interferometry in landslide complexes: insights from the Brienz/Brinzauls slope instability, Swiss Alps. Landslides, 2024 (https://doi.org/10.1007/s10346-024-02291-z)
- Louis, C., Halloran, L., Roques, C. Seasonal and Diurnal Freeze-Thaw Dynamics of a Rock Glacier and Their Impacts on Mixing and Solute Transport. EGUsphere preprint, 2024
- Singeisen, C., Massey, C., Wolter, A., Stahl, T. et al. Evolution of an earthquake-induced landslide complex in the South Island of New Zealand: How fault damage zones and seismicity contribute to slope failures. Geosphere, 2023 (https://doi.org/10.1130/GES02668.1)
- Hermle, D., Keuschnig, M., Krautblatter, M., Bickel, V.T. Systematic Quantification and Assessment of Digital Image Correlation Performance for Landslide Monitoring. Geosciences, 2023 (https://doi.org/10.3390/geosciences13120371)
- Morriss, M., Lehmann, B., Campforts, B., Brencher, G. et al. Alpine hillslope failure in the western US: insights from the Chaos Canyon landslide, Rocky Mountain National Park, USA. Earth Surface Dynamics (11), 2023 ([https://www.sciencedirect.com/science/article/pii/S0266352X23000253?via%3Dihub](https://esurf.copernicus.org/articles/11/1251/2023/))
- Hu, Y., Lu, Y. Study on soil-rock slope instability at mesoscopic scale using discrete element method. Computers and Geotechnics (157), 2023 (https://www.sciencedirect.com/science/article/pii/S0266352X23000253?via%3Dihub)
- Singeisen, C., Massey, C., Wolter, A., Kellett, R., Bloom, C., Stahl, T., Gasston, C., Jones, K. Mechanisms of rock slope failures triggered by the 2016 Mw 7.8 Kaikōura earthquake and implications for landslide susceptibility. Geomorphology (415), 2022 (https://doi.org/10.1016/j.geomorph.2022.108386)
- Bickel, V.T., Manconi, A. Decadal Surface Changes and Displacements in Switzerland. Journal of Geovisualization and Spatial Analysis 6(24), 2022 (https://doi.org/10.1007/s41651-022-00119-9)
- Xu, S., Fu, P., Quincey, D., Feng, M., Marsh, S., Lui, Q. UAV-based geomorphological evolution of the Terminus Area of the Hailuogou Glacier, Southeastern Tibetan Plateau between 2017 and 2020. Geomorphology (411), 2022 (https://doi.org/10.1016/j.geomorph.2022.108293)
- Aaron, J., Loew, S., Forrer, M. Recharge response and kinematics of an unusual earthflow in Liechtenstein. Landslides (18), 2021 (https://doi.org/10.1007/s10346-021-01633-5)
- Storni, E., Hugentobler, M., Manconi, A., Loew, S. Monitoring and analysis of active rockslide-glacier interactions (Moosfluh, Switzerland). Geomorphology (371), 2020 (https://doi.org/10.1016/j.geomorph.2020.107414)
- Raack, J., Conway, S.J., Heyer, T., Bickel, V.T., Philippe, P., Hiesinger, H., Johnsson, A., Massé, M. Present-day gully activity in Sisyphi Cavi, Mars – Flow-like features and block movements. Icarus (350), 2020 (https://doi.org/10.1016/j.icarus.2020.113899)
- Hermle, D., Keuschnig, M., Krautblatter, M. Potential of multisensor assessment using digital image correlation for landslide detection and monitoring. EGU 2020-16982, 2020 (https://doi.org/10.5194/egusphere-egu2020-16982)
- Storni, E., Loew, S., Manconi, A., Hugentobler, M. Monitoring and Analysis of Landslide-Glacier Interactions at the Great Aletsch Glacier (Switzerland). EGU 2020, 2020 (https://pdfs.semanticscholar.org/c735/e2d172970e3b8e83863a07a0c86a3bc3ec2a.pdf)
- Glueer, F., Loew, S., Manconi,, A., Aaron, J. From Toppling to Sliding: Progressive Evolution of the Moosfluh Landslide, Switzerland. JGR Earth Surface (124), 2019 (https://doi.org/10.1029/2019JF005019)
- Bickel, V.T., Manconi, A., Amann, F. Quantitative Assessment of Digital Image Correlation Methods to Detect and Monitor Surface Displacements of Large Slope Instabilities. Remote Sensing (10), 2018 (https://doi.org/10.3390/rs10060865)

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

- saveascii.m, by Carlos Adrian Vargas Aguilera.
Online: https://de.mathworks.com/matlabcentral/fileexchange/10243-saveascii-m-v6-2-oct-2009 (Status: May 22nd, 2020)

- gridfit.m, by John D'Errico.
Online: https://de.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit (Status: May 22nd, 2020)


### Hardware requirements:

- Input image size of 30 MB requires at least ~5 GB RAM


## Instructions

1. Clone and unzip DIC_FFT_ETHZ_Repo

2. Place primary and secondary image(s) in the "Input" folder (or copy images from 'Example' folder for a demo)

3. Customize setup for Preprocessing, DIC, and Postprocessing in run_pixel_offset.m, using the following variables. More information about the variables can be found further below:

          Inputs ----------------------------------------------------------------------------------
               - geotiff = 0/1
               - epsg = int.
               - coreg = 0/1
               - way = 1/2
               - yourtext = 'string'
               - pix = float.
               - master = 'string'
               - slave = 'string'

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

4. Execute run_pixel_offset.m

5. Collect displacement matrices ascii and deformation magnitude preview tif from "Output" folder for further utilization of results (e.g. drag & drop into a GIS). See "Graphical Abstract" for the organization of the output.

## Input Parameter Description
### Inputs   - - - - - - - - - -
#### geotiff
indicate here, if you would like to forward geotiff information (=1) or not (=0)

#### epsg
in case you would like to forward geotiff information, enter the desired EPSG code here
--> you can get this number here: https://spatialreference.org/ref/epsg/

#### coreg
indicate here, if you would like to manually select a region for the alignment (=1) or not (=0)

#### way
specify is you want the geotiff preview to be 8bit or 16bit

#### yourtext
insert a name for your output files

#### pix
insert the spatial resolution of the images used (GSD) in m/pixel

#### master (primary image)
Enter the name of the old image here, including the file type

#### slave (secondary image)
Enter the name of the new image here, including the file type

### Preprocessing (Wallis filter & Co-registration)   - - - - - - - - - -
#### wallis
indicate here, if you would like to apply a Wallis filter (=1) or not (=0)

#### win
specify Wallis filter window size, use even numbers only! A larger window will result in a larger context being considered for the dynamic contrast adaptation

#### tarm
specify the target mean value within each individual window - this mainly controls the image brightness (use the histogram of the filtered image to adapt the filter settings)

#### tars
specify the target standard deviation within each individual window - this mainly controls the radiometric normal distribution and contrast (use the histogram of the filtered image to adapt the filter settings)

#### b
specify the brightness enforcing constant (use the histogram of the filtered image to adapt the filter settings)

#### c
specify the contrast enforcing constant (use the histogram of the filtered image to adapt the filter settings)

#### sp
define the desired image split that is being used for the Co-Registration - 1  means that the entire image is matched, e.g. 8 means that the image is divided in 8 sub-parts, while every individual sub-image is being matched, where the mean correction of all 8 tiles will be used for Co-Registration

#### co_os
apply a oversampling factor for the Co-Registration, if desired - theoretically, a larger oversampling factor will result in a more accurate (sub-pixel) Co-Registration - in practice, a factor of 1 or max. 4 will be sufficient

### DIC   - - - - - - - - - -
#### wi
specify the desired window size in pixels; stick to the 'power of 2 rule' - a larger window will produce a better correlation and a better coverage of the displacement; a smaller window size helps to resolve smaller spatial scales of the displacement to the cost of a higher noise level; accuracy of the derived displacement is not influenced by the window size, see: http://www.mdpi.com/2072-4292/10/6/865
               
#### os
specify the oversampling factor of the FFT correlator - theoretically, a larger oversampling factor will result in a more accurate correlation (sub-pixel) - in practice, a factor of 1 or max. 4 will be sufficient, see http://www.mdpi.com/2072-4292/10/6/865

#### pix
enter the spatial resolution (GSD) of the used images in m/pixel

### Postprocessing (RMSE, Mean, Vector, Median filter)   - - - - - - - - - -
#### filter
choose whether you would like to use the RMSE filter (=1), the arithmetic mean filter (=2), the spatial vector filter (=3), or the median filter (=4) - the filters are described here: http://www.mdpi.com/2072-4292/10/6/865

#### thr
specify the RMSE threshold - a lower value will result in a strict filter, whereas a higher value will let pass values with a lower correlation quality (RMSE); a thr = 1 will result in a leaky filter, i.e., no filtering is performed at all - use this when you want to avoid any influence of a post-processing routine!

#### mfws
specify the mean filter window size in pixels - a larger window will result in a strongly blurred displacement map

#### cut
specify the mean filter cutoff value in pixels - this can be used to eliminate extreme noise that otherwise would mess with the mean filter
 
#### magcap
specify the magnitude (resultant) cap in pixels - values which are greater as this value are cut

#### xcap
specify the x direction (E-W, left-right) cap in pixels - values which are greater as this value are cut

#### ycap
specify the y direction (N-S, up-down) cap in pixels - values which are greater as this value are cut

#### med
specify the median filter window size in pixels - a larger window will result in a strongly blurred displacement map

### Additional settings   - - - - - - - - - -
#### scalax
specify the Min and Max values in x direction (E-W, left-right) that will be displayed in meters - everything above and below this range will be saturated at the given Min / Max value; in general, this range should be small for slow displacements and large for fast displacements; at the same time, regions with heterogeneous displacement velocities could require to this code with different ranges in order to catch the entire variety of displacement velocities!

#### scalay
specify the Min and Max values in y direction (N-S, up-down) that will be displayed in meters - everything above and below this range will be saturated at the given Min / Max value; in general, this range should be small for slow displacements and large for fast displacements; at the same time, regions with heterogeneous displacement velocities could require to this code with different ranges in order to catch the entire variety of displacement velocities!

#### coppia
provide a name for the ascii file that contains the displacement information (Output folder)

#### skip_x
specify the search window skip in x direction (E-W, left-right) - be aware of aliasing effects when going below the Nyquist limit (default)!

#### skip_y
specify the search window skip in y direction (N-S, up-down) - be aware of aliasing effects when going below the Nyquist limit (default)!

### Selection of results

<img src="https://github.com/bickelmps/DIC_FFT_ETHZ/blob/master/Figures/example2.png?raw=true">

##### Figure. (a) Optical image: Surface velocity at Moosfluh, Switzerland, measured with ground-based remote sensing methods. Average velocity obtained applying DIC on time lapse camera data (image pair 26–28 September 2018) acquired from across the valley. Unreliable and/or outlier DIC values are masked out. The large deformation occurring in the lower portions of the Moosfluh slope is well mapped from both techniques in a short temporal frame, and average velocities are in agreement with independent ground measurements. (b) Ground-based Radar Interferometry measurements (stack of 1 min interferograms over 8 h) performed at Moosfluh on 28 September 2017. Results shows good agreement with the results produced using DIC_FFT (a). Modified from Manconi et al., 2018.

<img src="https://github.com/bickelmps/DIC_FFT_ETHZ/blob/master/Figures/example1.png?raw=true">

##### Figure. DIC_FFT applied on Sentinel-1 and Sentinel-2 datasets over the Moosfluh unstable slope. White boundary delineates the approximate area of the Aletsch glacier and Moosfluh slope instability. (a) SAR intensity image: Resultant of the surface velocities measured across-track & along-track directions of Sentinel-1 in the period 8 June–18 October 2017. Background map is SLC image 8 June 2017 (slant range oriented towards north); (b) Multispectral image: Resultant of the horizontal surface velocities (E-N directions) measured on Sentinel-2 track R108 and R065 in the periods 19 July–17 October and 26 June–14 October 2017, respectively. Background maps are the surface reflectance for the images 19 July 2017 and 26 June 2017, respectively. Surface velocities in the lower portions of the Moosfluh instability and in the central part of the Great Aletsch glacier are well mapped by DIC. Inaccurate/outlier DIC values are masked out. Modified from Manconi et al., 2018.

------------------
##### MIT License - Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
