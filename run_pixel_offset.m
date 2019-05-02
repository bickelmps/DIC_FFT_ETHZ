%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%%_________________________________________________________________________
% for 8bit, equal dim, single- and three-band images, with on-demand
% geotiff information forwarding
% Pre-processing: Wallis filter, Co-Registration
% Post-processing: RMSE threshold, mean, median, spatial vector filter

% V. Bickel & A. Manconi 2.5.18
% valentin.bickel@erdw.ethz.ch / andrea.manconi@erdw.ethz.ch
% ETH Zurich / MPS Goettingen

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi

% Please cite this routine as:
% Bickel, V.T.; Manconi, A.; Amann, F.
% "Quantitative assessment of Digital Image Correlation methods to detect
% and monitor surface displacements of large slope instabilities."
% Remote Sens. 2018, 10(6), 865.
%%_________________________________________________________________________
%%
        warning off;clc
        clear
        clc
%%_________________________________________________________________________
%% INPUTS              
% Master image
geotiff = 1; % choose if you want to forward geotiff information, no = 0, yes = 1;
epsg = 21781; % EPSG code, specify in case you would like to forward geotiff information
PCS = 'CH1903 / LV03'; % PCS code, specify in case you would like to forward geotiff information
% Check here: https://spatialreference.org/ref/epsg/

if geotiff == 0
    cd Input
    master = imread('1_1.png'); % used for DIC
    orig_m = imread('1_1.png'); % used for plotting
    % Slave image
    slave = imread('1_6.png');
    orig_s = imread('1_6.png');
    cd ..
end
if geotiff == 1
    cd Input
    inputfilename = 'geotiff_2004';
    outfilename = 'geotiff_DIC_2004-2007';
    master = imread('geotiff_2004.tif'); % used for DIC
    orig_m = imread('geotiff_2004.tif'); % used for plotting
    % Slave image
    slave = imread('geotiff_2007.tif');
    orig_s = imread('geotiff_2007.tif');
    cd ..
end     
%%_________________________________________________________________________
%% Pre-processing: Wallis filter
cd DIC
disp '<<< Running! >>>'
wallis = 0; % wallis = 0, wallis filtering off, wallis = 1, wallis filtering on

win = 20; % block size = EVEN NUMBERS ONLY
tarm = 150; % target mean
tars = 150; % target standard deviation
b=1; % brightness enforcing constant
c=0.9995; % contrast enforcing constant

        if wallis == 1
        disp '<<< Pre-processing: Wallis filter >>>'
        tic
        [master,slave] = wallis_filter(master,slave,win,tarm,tars,b,c);
        toc
        end
        
clear tarm tars b c win
%%_________________________________________________________________________
%% Co-registration parameters
sp = 1;    % Image split for Co-registration
co_os = 1; % Image oversampling for Co-registration

%% DIC parameters            
% Offset type
wi = 64;   % Window size [pix]
os = 1;     % Oversampling factor
pix = 0.25; % Pixel size (in meters)

% Filter type
filter = 1; % 1 = threshold filter
            % 2 = arithmetic mean filter
            % 3 = vector filter
            % 4 = median filter
% #1
thr = 1;  % masking threshold (0-1, 1=no mask)

% #2
mfws = [6 6]; % mean filter window size for arithmetic mean filter
cut = wi; % window size/cut = cut off filter for a_mean filter, OPTIONAL

% #3
magcap = wi; % tolerance-diff, values which are greater as this value are cut [pixel]
xcap = wi/4; % tolerance-diff, values which are greater as this value are cut [pixel]
ycap = wi/4; % tolerance-diff, values which are greater as this value are cut [pixel]

% #4
med = [5 5]; % dimensions of median filter window [pixel]

% Colorscale -------------------------------------------
% min and max values for the X displacement colorscale, in meters
scalax=[-4 4];
% min and max values for the Y displacement colorscale, in meters
scalay=[-4 4];

% Additional definitions --------------------------------------------------
coppia = 't1-t0';    

% Additional PO analysis parameters
skip_x = wi/2; % half size of the window
skip_y = wi/2; % half size of the window

%%_________________________________________________________________________
%% Co-Registration
        disp '<<< Co-registration >>>'
        tic
        [T0,T1] = co_registration(master,slave,sp,co_os);        
        toc         
%%_________________________________________________________________________
%% DIC
        disp '<<< Running FFT >>>'
        tic
        [R] = pixoff(T0,T1,skip_x,skip_y,wi,os,coppia);

        cd ..
        cd Output
        %t=load(['pr_',coppia,'.txt']); % load pixoff results [pix]
        t = R;
        cd .. 
        cd DIC
        ii=1:numel(min(t(:,1)):skip_x:max(t(:,1))); jj=1:numel(min(t(:,2)):skip_x:max(t(:,2)));
        [I,J]=meshgrid(ii,jj);
        i=reshape(I,[numel(I),1]); j=reshape(J,[numel(J),1]);
        nx = round(size(T0,2)/5); ny = round(size(T0,1)/10);
        toc       
%%_________________________________________________________________________
%% Filtering
        disp '<<< Post-processing: Filtering >>>'
        tic
        %_#1_______________________________________________________________
        % Threshold correlation filter
        if filter == 1
        filter_selection = 1;
        pp = t(:,end)<thr; % gives a 1 true or 0 not true
        end
        %__________________________________________________________________
        %_#2_______________________________________________________________
        % Arithmetic mean filter
        if filter == 2
        filter_selection = 2;
        [tpostfiltx,tpostfilty] = mean_filt(T0,t,mfws,wi);
        end
        %__________________________________________________________________
        %_#3_______________________________________________________________
        % Vector filter
        if filter == 3
        filter_selection = 3;
        [tot_check] = vec_filt(t,T0,wi,ycap,xcap,magcap);
        end
        %__________________________________________________________________
        %_4________________________________________________________________
        % Median filter
        if filter == 4
        filter_selection = 2;
        [tpostfiltx,tpostfilty] = med_filt(T0,t,med,wi);
        end
        %__________________________________________________________________
        toc
%%_________________________________________________________________________
%% Forwarding Geotiff information
% Credits to F. Gluer & M. Haeusler, SED Zurich
if geotiff == 1
    [a,b] = size(J);
    outfigure = reshape((sqrt(t(pp,3).^2+t(pp,4).^2)*pix),a,b);
    outfigure=outfigure';
    outfigure(:,end)=outfigure(:,a-1);

    [RR] = geotiffinfo(sprintf('../Input/%s.tif',inputfilename));
    %Create worldfile .tfw
    P1=[RR.CornerCoords.X(4) RR.CornerCoords.Y(4)]; % P1: coordinates lower-left corner
    P2=[RR.CornerCoords.X(3) RR.CornerCoords.Y(3)]; % Pn: coordinates [...]         
    P4=[RR.CornerCoords.X(1) RR.CornerCoords.Y(1)];

    r1 = sqrt((P1(1)-P4(1))^2+(P1(2)-P4(2))^2);                        
    beta2 = 2*asin(sqrt((P1(1)-P4(1))^2+((P1(2)+r1)-P4(2))^2)/(2*r1)); % rotation angle of image
    r2 = sqrt((P1(1)-P2(1))^2+(P1(2)-P2(2))^2);
    alpha2 = 2*asin(sqrt(((P1(1)+r2)-P2(1))^2+(P1(2)-P2(2))^2)/(2*r2)); % rotation angle of image

    size_pixel_x = r2/a; % distance of one pixel in x-direction [m]
    size_pixel_y = r1/b; % distance of one pixel in y-direction [m]

    % https://en.wikipedia.org/wiki/World_file
    %     Line 1: A: x-component of the pixel width (x-scale)
    %     Line 2: D: y-component of the pixel width (y-skew)
    %     Line 3: B: x-component of the pixel height (x-skew)
    %     Line 4: E: y-component of the pixel height (y-scale), typically negative
    %     Line 5: C: x-coordinate of the center of the original image's upper left pixel transformed to the map
    %     Line 6: F: y-coordinate of the center of the original image's upper left pixel transformed to the map

    fid = fopen(sprintf('../Output/%s_worldfile.tfw',inputfilename), 'wt'); % write to .tfw worldfile
    fprintf(fid, '%10.10f\n', size_pixel_x*cos(alpha2));
    fprintf(fid, '%10.10f\n', -size_pixel_x*sin(alpha2));
    fprintf(fid, '%10.10f\n', -size_pixel_y*sin(beta2));
    fprintf(fid, '%10.10f\n', -size_pixel_y*cos(beta2));
    fprintf(fid, '%10.10f\n', P4(1));
    fprintf(fid, '%10.10f\n', P4(2));
    fclose(fid);

    R2 = worldfileread(sprintf('../Output/%s_worldfile.tfw',inputfilename),'planar',[b, a]);
    if RR.PCS == 'CH1903 / LV03' % Swiss coordinate system
        geotiffwrite(sprintf('../Output/%s',outfilename),outfigure,R2,'CoordRefSysCode',epsg);
    else
        try
            RR.PCS == PCS; % desired coordinate system
            geotiffwrite(sprintf('../Output/%s',outfilename),outfigure,R2,'CoordRefSysCode',epsg);
        catch
            print('Coordinate system not defined - please specify EPSG-code and PCS')
        end
    end
end
%%_________________________________________________________________________
%% Plotting
% Multi-band reduction of input image, if required
        [x, y, z] = size(orig_m);
        if z > 3
            orig_m_int(:,:,1) = orig_m(:,:,1);
            orig_m_int(:,:,2) = orig_m(:,:,2);
            orig_m_int(:,:,3) = orig_m(:,:,3);
            orig_m = orig_m_int;
        
            orig_s_int(:,:,1) = orig_s(:,:,1);
            orig_s_int(:,:,2) = orig_s(:,:,2);
            orig_s_int(:,:,3) = orig_s(:,:,3);
            orig_s = orig_s_int;
        end
% Plot 2D offset and displacement vectors 
        if filter_selection == 1 % RMSE threshold filter __________________
        figure(1)
        
        h = subplot(2,2,1);
        imshow(orig_m)
        caxis([-100 180])
        colormap(h,bone);
        title('Master image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
                
        i = subplot(2,2,2);
        imshow(orig_s)
        caxis([-100 180])
        colormap(i,bone);
        title('Slave image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
               
        g = subplot(2,2,3);
        twoD=sqrt(t(pp,3).^2+t(pp,4).^2);
        scatter(t(pp,2),t(pp,1),30,(twoD.*pix),'filled','s','MarkerFaceAlpha',.95,'MarkerEdgeAlpha',.95)
        colormap(jet);
        title('2D Offset magnitude')
        caxis(scalay)
        xlim([0 size(T0,2)])
        ylim([0 size(T0,1)])
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
        pbaspect([size(T0,2) size(T0,1) 1])
        set(gca,'YDir','reverse')
        gc=colorbar('EastOutside');
        xlabel(gc,'meters')
        
        m = subplot(2,2,4);
        quiver(t(:,2),t(:,1),(t(:,4)*pix),(t(:,3)*pix),5,'b')
        colormap(jet);
        title('Displacement vectors')
        xlim([0 size(T0,2)])
        ylim([0 size(T0,1)])
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
        pbaspect([size(T0,2) size(T0,1) 1])
        set(gca,'YDir','reverse')
        
        linkaxes([h,g,i,m], 'xy')
        end
        if filter_selection == 2 % Mean & Median filter ___________________
        figure(1)
        
        h = subplot(2,2,1);
        imshow(orig_m)
        caxis([-100 180])
        colormap(h,bone);
        title('Master image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
                
        i = subplot(2,2,2);
        imshow(orig_s)
        caxis([-100 180])
        colormap(i,bone);
        title('Slave image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
               
        g = subplot(2,2,3);
        twoD=sqrt(tpostfiltx.^2+tpostfilty.^2);
        scatter(t(:,2),t(:,1),30,(twoD*pix),'filled','s','MarkerFaceAlpha',.95,'MarkerEdgeAlpha',.95)       
        colormap(jet);
        title('2D Offset magnitude')
        caxis(scalay)
        xlim([0 size(T0,2)])
        ylim([0 size(T0,1)])
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
        pbaspect([size(T0,2) size(T0,1) 1])
        set(gca,'YDir','reverse')
        gc=colorbar('EastOutside');
        xlabel(gc,'meters')
        
        m = subplot(2,2,4);
        quiver(t(:,2),t(:,1),(tpostfiltx(:,1)*pix),(tpostfilty(:,1)*pix),5,'b')
        colormap(jet);
        title('Displacement vectors')
        xlim([0 size(T0,2)])
        ylim([0 size(T0,1)])
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
        pbaspect([size(T0,2) size(T0,1) 1])
        set(gca,'YDir','reverse')
        
        linkaxes([h,g,i,m], 'xy')     
        end
        if filter_selection == 3 % Vector filter __________________________
        figure(1)
        
        h = subplot(2,2,1);
        imshow(orig_m)
        caxis([-100 180])
        colormap(h,bone);
        title('Master image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
                
        i = subplot(2,2,2);
        imshow(orig_s)
        caxis([-100 180])
        colormap(i,bone);
        title('Slave image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
               
        g = subplot(2,2,3);
        scatter(t(:,2),t(:,1),30,(tot_check(:,5)*pix),'filled','s','MarkerFaceAlpha',.95,'MarkerEdgeAlpha',.95)
        colormap(jet);
        title('2D Offset magnitude')
        caxis(scalay)
        xlim([0 size(T0,2)])
        ylim([0 size(T0,1)])
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
        pbaspect([size(T0,2) size(T0,1) 1])
        set(gca,'YDir','reverse')
        gc=colorbar('EastOutside');
        xlabel(gc,'meters')
        
        m = subplot(2,2,4);
        quiver(t(:,2),t(:,1),(tot_check(:,4)*pix),(tot_check(:,3)*pix),5,'b')
        colormap(jet);
        title('Displacement vectors')
        xlim([0 size(T0,2)])
        ylim([0 size(T0,1)])
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
        pbaspect([size(T0,2) size(T0,1) 1])
        set(gca,'YDir','reverse')
        
        linkaxes([h,g,i,m], 'xy')     
        end
        cd ..
%% END OF SCRIPT        
% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi