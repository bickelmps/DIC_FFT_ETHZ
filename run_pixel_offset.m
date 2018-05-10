%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%%_________________________________________________________________________
% for 8bit, equal dim, single band images
% Pre-processing: Wallis filter
% Post-processing: RSME threshold, mean, median, vector filter

% A. Manconi & V. Bickel, 2.5.18
% andrea.manconi@erdw.ethz.ch / valentin.bickel@erdw.ethz.ch
% ETH Zurich / MPS Goettingen

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi

% Please cite this routine as:
% Bickel, V.T.; Manconi, A.; Amann, F.
% "Quantitative assessment of Digital Image Correlation methods to detect
% and monitor surface displacements of large slope instabilities."
% Remote Sensing Journal 2018.
%%_________________________________________________________________________
%%
        warning off;clc
        clear
        clc
%%_________________________________________________________________________
%% INPUTS              
% Master image
cd Input
master = imread('INPUT_IMAGE_1.type'); % used for DIC
orig_m = imread('INPUT_IMAGE_1.type'); % used for plotting
% Slave image
slave = imread('INPUT_IMAGE_2.type');
orig_s = imread('INPUT_IMAGE_2.type');
cd ..
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
wi = 128;   % Window size [pix]
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
xcap = wi/8; % tolerance-diff, values which are greater as this value are cut [pixel]
ycap = wi/8; % tolerance-diff, values which are greater as this value are cut [pixel]

% #4
med = [5 5]; % dimensions of median filter window [pixel]

% Colorscale -------------------------------------------
% min and max values for the X displacement colorscale, in meters
scalax=[-5 5];
% min and max values for the Y displacement colorscale, in meters
scalay=[-5 5];

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
        t=load(['pr_',coppia,'.txt']); % load pixoff results [pix]
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
%% Plotting
        % Plot 2D offset and displacement vectors 
        if filter_selection == 1 % RSME threshold filter __________________
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
        
% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi