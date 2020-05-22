function [tot_check] = vec_filt(t,T0,wi,ycap,xcap,magcap)
%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%% GIT VERSION
%% VECTOR FILTER
%%_________________________________________________________________________
%{
for 8bit, equal dim, single- and three-, and four-band images, with
GIS-ready geotiff preview

*Pre-processing*: Wallis filter, Co-Registration
*Post-processing*: RMSE threshold, mean, median, spatial vector filter
___________________________________________________________________________
STRUCTURE:
|PARAMETERS - modify
         |CODE - don't modify
___________________________________________________________________________
V. Bickel & A. Manconi 21.5.2020
valentin.bickel@erdw.ethz.ch / andrea.manconi@erdw.ethz.ch
ETH Zurich / MPS Goettingen
---------------------------------------------------------------------------
MIT License
Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
---------------------------------------------------------------------------
Please cite this routine as:
Bickel, V.T.; Manconi, A.; Amann, F.
"Quantitative assessment of Digital Image Correlation methods to detect
and monitor surface displacements of large slope instabilities."
Remote Sens. 2018, 10(6), 865.
%}
%%_________________________________________________________________________
%%
% t = output pixoff.m or pixoff_ncc.m (deformation matrix nx5)
% T0 = input master image (e.g. mxn uint8)
% wi = pixoff search window size (num)
% ycap = max allowed difference between mean y (3x3)and local y deformation
% xcap = max allowed difference between mean x (3x3) and local x deformation
% magcap = max allowed difference between mean mag (3x3) and local mag deformation

% OUTPUT
% tot_check = vector filtered deformation matrix t (nx5)

        % Preparation
        t_res = sqrt((t(:,3).^2) + (t(:,4).^2));% Resultants of t
        size_tx = ceil(size(T0)/(wi/2));
        size_ty = size(t,1);
        hsize = [3 3]; % mean filter size (3x3 was determined as optimum)
        h = fspecial('average', hsize);
        t1 = t(:,3); % y values
        Mt1 = reshape(t1,size_tx(1:2)); % 2D Matrix of y def
        % Mt1(Mt1>(wi/cut))=0;      % OPTIONAL CUT OFF FILTER FOR OUTLIERS 
        % Mt1(Mt1<(wi/cut*(-1)))=0; % OPTIONAL CUT OFF FILTER FOR OUTLIERS
        t2 = t(:,4); % x values
        Mt2 = reshape(t2,size_tx(1:2)); % 2D Matrix of x def
        % Mt2(Mt2>(wi/cut))=0;      % OPTIONAL CUT OFF FILTER FOR OUTLIERS 
        % Mt2(Mt2<(wi/cut*(-1)))=0; % OPTIONAL CUT OFF FILTER FOR OUTLIERS
        t3 = t_res(:,1); % Magnitude values
        Mt3 = reshape(t3,size_tx(1:2)); % 2D Matrix of resultants
        
        % Filtering
        % X values
        Imf1 = imfilter(Mt1,h); % filter
        Imf1_T = Imf1'; % Reshaping
        imf1 = Imf1_T(:);
        imf1 = imf1(1:size_ty); % cutted for plotting
        abs1 = abs(t1-imf1); % preparation of comparison vector
        abs1(abs1>ycap) = NaN; % threshold cap filtering
        compare_y = cat(2,t1,abs1);
        compare_y(any(isnan(compare_y),2),:)=NaN;
        abs1out = compare_y(:,1);
        
        % Y values
        Imf2 = imfilter(Mt2,h); % filter
        Imf2_T = Imf2'; % Reshaping
        imf2 = Imf2_T(:);
        imf2 = imf2(1:size_ty); % cutted for plotting
        abs2 = abs(t2-imf2); % preparation of comparison vector
        abs2(abs2>xcap) = NaN; % threshold cap filtering
        compare_x = cat(2,t2,abs2);
        compare_x(any(isnan(compare_x),2),:)=NaN;
        abs2out = compare_x(:,1);
        
        % Magnitude values
        Imf3 = imfilter(Mt3,h); % filter
        Imf3_T = Imf3'; % Reshaping
        imf3 = Imf3_T(:);
        imf3 = imf3(1:size_ty); % cutted for plotting
        abs3 = abs(t3-imf3); % preparation of comparison vector
        abs3(abs3>magcap) = NaN; % threshold cap filtering
        compare_mag = cat(2,t3,abs3);
        compare_mag(any(isnan(compare_mag),2),:)=NaN;
        abs3out = compare_mag(:,1);
        
        % Finalization & Output 
        tot_check = cat(2,t(:,1:2),abs1out,abs2out,abs3out); % use magnitude as additional filter
        tot_check(any(isnan(tot_check),2),:)=0; % eliminate every pixel with modified data
        
        cd Output
        save(['pr_t1-t0','.txt'],'tot_check', '-ascii');
        cd ..
%%
%{
MIT License
Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
---------------------------------------------------------------------------
Please cite this routine as:
Bickel, V.T.; Manconi, A.; Amann, F.
"Quantitative assessment of Digital Image Correlation methods to detect
and monitor surface displacements of large slope instabilities."
Remote Sens. 2018, 10(6), 865.
%}     