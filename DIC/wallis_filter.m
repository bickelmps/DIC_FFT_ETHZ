function [master,slave] = wallis_filter(master,slave,win,tarm,tars,b,c)
%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%% GIT VERSION
%% WALLIS FILTER
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
img1 = im2double(master);
img2 = im2double(slave);
bit = 255;

% Padding
im1 = padarray(img1,[win/2 win/2],'symmetric');
im2 = padarray(img2,[win/2 win/2],'symmetric');

% Loop img dimensions
lenx = size(im1,1);
leny = size(im1,2);

% imfilter of im
hsize = [win+1 win+1];
h = fspecial('average', hsize);
imf1 = imfilter(im1,h);
imf2 = imfilter(im2,h);

% stdfilt of im
imstd1 = stdfilt(im1,ones(21));
imstd2 = stdfilt(im2,ones(21));

for i = 1:lenx
    for j = 1:leny;
        imG1(i,j)= (im1(i,j)-imf1(i,j))*c*tars/(c*imstd1(i,j)+(1-c)*tars)+b*tarm+(1-b)*imf1(i,j);
        if imG1(i,j)<0
            imG1(i,j)=0;
        elseif imG1(i,j)>255
            imG1(i,j)=255;
        end
    end 
    disp(i)
end
for i = 1:lenx
    for j = 1:leny;
        imG2(i,j)= (im2(i,j)-imf2(i,j))*c*tars/(c*imstd2(i,j)+(1-c)*tars)+b*tarm+(1-b)*imf2(i,j);
        if imG2(i,j)<0
            imG2(i,j)=0;
        elseif imG2(i,j)>255
            imG2(i,j)=255;
        end
    end  
    disp(i)
end

% De-padding
imG1 = imG1(1+win/2:end-win/2,1+win/2:end-win/2);
master = uint8(imG1);
imG2 = imG2(1+win/2:end-win/2,1+win/2:end-win/2);
slave = uint8(imG2);
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