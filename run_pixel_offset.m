%%_________________________________________________________________________
%% DIC using FFT-approach and optional pre- and post-processing
%% GIT VERSION
%% MAIN
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
        warning off;clc
        clear
        clc
%%_________________________________________________________________________
%% INPUTS              
tic
coreg = 0; % 0 = full image, 1 = define patch
geotiff = 0; % choose if you want to output a geotiff preview, no = 0, yes = 1;
way = 1; % define way to output quickview getiff, default 8bit output = 1 , alternative 16bit output = 2. NOTE: choice changes unit of preview!!!
epsg = 21781; % EPSG code, specify in case you would like to output geotiff preview
yourtext = 'dic_run_00'; % specify the name of the output
pix = 0.25; % Pixel size or GSD (in meters), only required if geotiff = 0
  
        cd Input
        if geotiff == 1
        % Primary image
[primary, R] = geotiffread('YOURIMAGEHERE.tif'); % read
        % Secondary image
[secondary, ~]  = geotiffread('YOURIMAGEHERE.tif'); % read
         pix = R.CellExtentInWorldY; % GSD in meters/pixel taken from the geotiff
         latcoord = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
         loncoord = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
         
        else
        %Primary image
[primary] = imread('YOURIMAGEHERE.tif'); %  read
        %Secondary image
[secondary] = imread('YOURIMAGEHERE.tif'); % read
        end
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
        [primary,secondary] = wallis_filter(primary,secondary,win,tarm,tars,b,c);
        toc
        end
        
clear wallis tarm tars b c win
%%_________________________________________________________________________
%% Co-registration parameters
sp = 1;    % Image split for Co-registration
co_os = 1; % Image oversampling for Co-registration

%% DIC parameters            
% Offset type
wi = 128;   % Window size [pix]
os = 1;     % Oversampling factor

% Filter type
filter = 1; % 1 = threshold filter
            % 2 = arithmetic mean filter
            % 3 = vector filter
            % 4 = median filter
% #1
thr = 1;  % masking threshold (0-1, 1=no mask)

% #2
mfws = [6 6]; % mean filter window size for arithmetic mean filter
cut = 12; % window size/cut = cut off filter for a_mean filter, OPTIONAL

% #3
magcap = wi; % tolerance-diff, values which are greater as this value are cut [pixel]
xcap = wi/4; % tolerance-diff, values which are greater as this value are cut [pixel]
ycap = wi/4; % tolerance-diff, values which are greater as this value are cut [pixel]

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
skip_x = wi/2; % half size of the window = Nyquist is happy!
skip_y = wi/2;
skip = skip_x;
%%_________________________________________________________________________
%% Co-Registration
        if coreg == 1
         imagesc(primary(:,:,1)); colormap gray; title("NOTE: if bbox size is too small, coreg might fail")
         set(gcf,'position',get(0,'ScreenSize'));
         rect = getrect(gca);
         close
        else
         rect = 1;
        end
        
        disp '<<< Co-registration >>>'
        tic
        [T0,T1] = co_registration(primary,secondary,sp,co_os,coreg,rect);        
        toc  
        clear sp co_os coreg
%%_________________________________________________________________________
%% DIC
        disp '<<< Running FFT >>>'
        tic
        [RR] = pixoff(T0,T1,skip_x,skip_y,wi,os,coppia);

        clear os
        
        % txt ASCII output 
        ii = (min(RR(:,1)):skip:max(RR(:,1)));
        jj = (min(RR(:,2)):skip:max(RR(:,2)));

        DX = RR(:,4).*pix; % left right displacement (E W) + pixel to meter
        DY = RR(:,3).*pix; % up down displacement (N S) + pixel to meter
        D2D = sqrt(DX.^2+DY.^2); % displacement resultant (magnitude)
        RMSE = RR(:,5); % RMSE taken from the FFT correlation
        
        % orientation correction, because Matlab is weird
        DX_mat = flipud(rot90(reshape(DX, [size(jj,2),size(ii,2)]))); % reshape to matrix and correct
        DY_mat = flipud(rot90(reshape(DY, [size(jj,2),size(ii,2)])));
        D2D_mat = flipud(rot90(reshape(D2D, [size(jj,2),size(ii,2)])));
        RMSE_mat = flipud(rot90(reshape(RMSE, [size(jj,2),size(ii,2)])));
        
        DX = reshape(flipud(DX_mat)',size(D2D)); % reshape to vector
        DY = reshape(flipud(DY_mat)',size(D2D));
        D2D = reshape(flipud(D2D_mat)',size(D2D));
        RMSE = reshape(flipud(RMSE_mat)',size(D2D));
        
        if geotiff == 1
            res = [loncoord(RR(:,2))',latcoord(RR(:,1))',DX,DY,D2D,RMSE];
        else
            res = [RR(:,2),RR(:,1),DX,DY,D2D,RMSE];
        end
        cd ..
        addpath(genpath('DIC'))
        cd Output
        saveascii(res,[yourtext,'.txt'],2, ','); % LONcoord | LATcoord | xoff(m) | yoff(m) | horizontal resultant (m) | RMSE
        cd ..
        t = RR;
        nx = round(size(T0,2)/5); ny = round(size(T0,1)/10);
        
        clear ii jj DX DY RMSE D2D DX_mat DY_mat D2D_mat RMSE_mat res

        ii=1:numel(min(t(:,1)):skip_x:max(t(:,1))); jj=1:numel(min(t(:,2)):skip_x:max(t(:,2)));
        [I,J]=meshgrid(ii,jj);
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
        disp '<<< Post-processing: geotiff preview >>>'
        tic
        if geotiff == 1
         [a,b] = size(J);
         if filter_selection == 1
          pp(pp == 0) = 1; % vector correction for reshape
          outfigure = reshape((sqrt(t(pp,3).^2+t(pp,4).^2)*pix),a,b);
          outfigure=outfigure';
          outfigure(:,end)=outfigure(:,a-1);
         end
        if filter_selection == 2
         outfigure = reshape((sqrt(tpostfiltx.^2+tpostfilty.^2)*pix),a,b);
         outfigure=outfigure';
         outfigure(:,end)=outfigure(:,a-1);       
        end
        if filter_selection == 3
         outfigure = reshape((tot_check(:,5)*pix),a,b);
         outfigure=outfigure';
         outfigure(:,end)=outfigure(:,a-1);
        end
   
         load mycmap_s.mat % custom colormap, choose if custom or Matlab map is used for the preview
         %mycmap = colormap(jet); % Matlab colormap
         %close
         Q2 = imresize(outfigure,size(primary(:,:,1)),'near');

         cd Output
         if way == 1
         geotiffwrite([yourtext,'.tif'], uint8(Q2), mycmap,R,'CoordRefSysCode',epsg);
         end
         if way == 2 
         geotiffwrite([yourtext,'.tif'], uint16(Q2*1000),R,'CoordRefSysCode',epsg); % factor 1000 to maintain small signals in the 16bit output, adapt as required
         % NOTE: the original output is in m/pix, i.e., applying a factor of e.g. *1000  will result in a unit of mm/pix for the preview output !!!
         end
         cd ..     
        end
        toc
%%_________________________________________________________________________
%% Plotting
% Multi-band reduction of input images, if required
        [x, y, z] = size(primary);
        if z > 3
            orig_m_int(:,:,1) = primary(:,:,1);
            orig_m_int(:,:,2) = primary(:,:,2);
            orig_m_int(:,:,3) = primary(:,:,3);
            orig_m = orig_m_int;
        
            orig_s_int(:,:,1) = secondary(:,:,1);
            orig_s_int(:,:,2) = secondary(:,:,2);
            orig_s_int(:,:,3) = secondary(:,:,3);
            orig_s = orig_s_int;
        else
            orig_m = primary;
            orig_s = secondary;
        end
        clear primary secondary orig_m_int orig_s_int
% Plot 2D offset and displacement vectors 
        if filter_selection == 1 % RMSE threshold filter __________________
        figure(1)
        
        h = subplot(2,2,1);
        imshow(orig_m)
        caxis([-100 180])
        colormap(h,gray);
        title('primary image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
                
        i = subplot(2,2,2);
        imshow(orig_s)
        caxis([-100 180])
        colormap(i,gray);
        title('secondary image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
               
        g = subplot(2,2,3);
        twoD=sqrt(t(pp,3).^2+t(pp,4).^2);
        scatter(t(pp,2),t(pp,1),30,(twoD.*pix),'filled','s','MarkerFaceAlpha',.95,'MarkerEdgeAlpha',.95)
        colormap(g,jet);
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
        xlabel(gc,'pixel units')
        
        m = subplot(2,2,4);
        quiver(t(:,2),t(:,1),(t(:,4)*pix),(t(:,3)*pix),5,'b')
        colormap(m,jet);
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
        colormap(h,gray);
        title('primary image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
                
        i = subplot(2,2,2);
        imshow(orig_s)
        caxis([-100 180])
        colormap(i,gray);
        title('secondary image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
               
        g = subplot(2,2,3);
        twoD=sqrt(tpostfiltx.^2+tpostfilty.^2);
        scatter(t(:,2),t(:,1),30,(twoD*pix),'filled','s','MarkerFaceAlpha',.95,'MarkerEdgeAlpha',.95)       
        colormap(g,jet);
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
        xlabel(gc,'pixel units')
        
        m = subplot(2,2,4);
        quiver(t(:,2),t(:,1),(tpostfiltx(:,1)*pix),(tpostfilty(:,1)*pix),5,'b')
        colormap(m,jet);
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
        colormap(h,gray);
        title('primary image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
                
        i = subplot(2,2,2);
        imshow(orig_s)
        caxis([-100 180])
        colormap(i,gray);
        title('secondary image')
        box on; axis on; grid on; grid minor
        ax = gca;
        ax.GridAlpha=0.8;
        xticks([nx 2*nx 3*nx 4*nx 5*nx])
        yticks([ny 2*ny 3*ny 4*ny 5*ny 6*ny 7*ny 8*ny 9*ny 10*ny])
               
        g = subplot(2,2,3);
        scatter(t(:,2),t(:,1),30,(tot_check(:,5)*pix),'filled','s','MarkerFaceAlpha',.95,'MarkerEdgeAlpha',.95)
        colormap(g,jet);
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
        xlabel(gc,'pixel units')
        
        m = subplot(2,2,4);
        quiver(t(:,2),t(:,1),(tot_check(:,4)*pix),(tot_check(:,3)*pix),5,'b')
        colormap(m,jet);
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
        
%% END OF SCRIPT        
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