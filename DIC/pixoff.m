function [R]=pixoff(An,Bn,skip_x,skip_y,wi,os,coppia)
%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%% GIT VERSION
%% RUN PIXEL OFFSET
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
An=padarray(An,[wi/2, wi/2]);
Bn=padarray(Bn,[wi/2, wi/2]);
sA=size(An);     
    start_x=1;
    end_x=sA(2);
    start_y=1;
    end_y=sA(1);
k=1;
for i=start_y+((wi/2)):skip_y:end_y-((wi/2))    
    for j=start_x+((wi/2)):skip_x:end_x-((wi/2))        
% Chips
P1=An(i-(wi/2):i+(wi/2),j-(wi/2):j+(wi/2));
P2=Bn(i-(wi/2):i+(wi/2),j-(wi/2):j+(wi/2));
    output = dftregistration(fft2(P1),fft2(P2),os);
%%%Result in coordinates
R(k,:)=[i-wi/2 j-wi/2 -output(3) -output(4) output(1)];
    k=k+1;
    end      
end
% Write results
cd ..
cd Output
     save(['pr_',coppia,'.txt'],'R', '-ascii'); % x position | y position | xoff(pix) | yoff(pix) | RMSE
cd ..
cd DIC
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