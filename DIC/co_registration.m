function [T0,T1] = co_registration(master,slave,sp,co_os)
%%_________________________________________________________________________
%% Image Co-registration
% A. Manconi & V. Bickel, 2.5.18
% andrea.manconi@erdw.ethz.ch / valentin.bickel@erdw.ethz.ch
% ETH Zurich / MPS Goettingen

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi
%%_________________________________________________________________________
%%
% convert to singleband, if not already single band
% T0=rgb2gray(A);
% T1=rgb2gray(B);
T0=master;
T1=slave;
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
T1=shift(double(T1),round(mean(out(:,3)))  ,round(mean(out(:,4))));

% MIT License
% Copyright (c) 2018 Valentin Tertius Bickel & Andrea Manconi