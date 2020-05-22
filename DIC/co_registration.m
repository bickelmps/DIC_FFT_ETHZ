function [T0,T1] = co_registration(master,slave,sp,co_os,coreg,rect)
%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%% GIT VERSION
%% CO-REGISTRATION
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
% convert to singleband, if not already single band
% T0=rgb2gray(A);
% T1=rgb2gray(B);
if coreg == 1
    T0=imcrop(master(:,:,1:3), rect);
    T1=imcrop(slave(:,:,1:3), rect);
else
    T0=master;
    T1=slave;
end
% find indexes
pp_r=1:size(T0,1)/sp:size(T0,1);
pp_c=1:size(T0,2)/sp:size(T0,2);
pp_r=[round(pp_r) size(T0,1)];
pp_c=[round(pp_c) size(T0,2)];
kk=1;
for i=1:sp    
    for j=1:sp        
     % Window 
co_w_xs=pp_r(i);
co_w_xe=pp_r(i+1);
co_w_ys=pp_c(j);
co_w_ye=pp_c(j+1);   
[out(kk,:)]=dftregistration(fft2(T0(co_w_xs:co_w_xe,co_w_ys:co_w_ye)),...
                              fft2(T1(co_w_xs:co_w_xe,co_w_ys:co_w_ye)),co_os);                                    
                          kk=kk+1;
    end
end
T1=shift(double(slave),round(mean(out(:,3)))  ,round(mean(out(:,4))));
T0 = master;
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