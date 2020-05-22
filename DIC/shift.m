function [B]=shift(A,nr,nc)
%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%% GIT VERSION
%% SHIFT
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
[M,N]=size(A);
B=zeros(M,N);

if nr>=0 & nc>=0
    B(nr+1:M,nc+1:N)=A(1:M-nr,1:N-nc);
elseif nr<0 & nc<0
    B(1:M+nr,1:N+nc)=A(1-nr:M,1-nc:N);
elseif nr>=0 & nc<0
    B(nr+1:M,1:N+nc)=A(1:M-nr,1-nc:N);
elseif nr<0 & nc>=0
    B(1:M+nr,nc+1:N)=A(1-nr:M,1:N-nc); 
end
    
return
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