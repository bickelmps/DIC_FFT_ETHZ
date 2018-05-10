function [master,slave] = wallis_filter(master,slave,win,tarm,tars,b,c)
%%_________________________________________________________________________
%% Wallis filter
% A. Manconi & V. Bickel, 2.5.18
% andrea.manconi@erdw.ethz.ch / valentin.bickel@erdw.ethz.ch
% ETH Zurich / MPS Goettingen

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi

% Underlying principle taken from:
% Baltsavias, E. P. "Multiphoto geometrically constrained matching."
% Doctoral Thesis, ETH Zurich, (1991).
%%_________________________________________________________________________
%% Load image
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

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi