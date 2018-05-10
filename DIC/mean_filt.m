function [tpostfiltx,tpostfilty] = mean_filt(T0,t,mfws,wi);
%%_________________________________________________________________________
%% Mean filter for 2D displacements
% A. Manconi & V. Bickel, 2.5.18
% andrea.manconi@erdw.ethz.ch / valentin.bickel@erdw.ethz.ch
% ETH Zurich / MPS Goettingen

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
%%_________________________________________________________________________
%%
        size_tx = size(T0,2); size_ty = size(t,1); Mtfilter = fspecial('average',mfws);
        t1 = t(:,3); Mt1 = vec2mat(t1,round((size_tx/(wi/2))));
        % Mt1(Mt1>(wi/cut))=0; Mt1(Mt1<(wi/cut*(-1)))=0; % OPTIONAL CUT OFF FILTER FOR OUTLIERS 
        Mtpostfilty = filter2(Mtfilter,Mt1); Mtpostfilty_T = Mtpostfilty'; tpostfilty = Mtpostfilty_T(:);
        t2 = t(:,4); Mt2 = vec2mat(t2,round((size_tx/(wi/2))));
        % Mt2(Mt2>(wi/cut))=0; Mt2(Mt2<(wi/cut*(-1)))=0; % OPTIONAL CUT OFF FILTER FOR OUTLIERS 
        Mtpostfiltx = filter2(Mtfilter,Mt2); Mtpostfiltx_T = Mtpostfiltx'; tpostfiltx = Mtpostfiltx_T(:);
        tpostfilty = tpostfilty(1:size_ty); tpostfiltx = tpostfiltx(1:size_ty);
        
        cd ..
        cd Output
        Rmf = cat(2,t(:,1),t(:,2),tpostfilty(:,1),tpostfiltx(:,1),t(:,5));
        save(['pr_t1-t0','.txt'],'Rmf', '-ascii');
        cd ..
        cd DIC

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi