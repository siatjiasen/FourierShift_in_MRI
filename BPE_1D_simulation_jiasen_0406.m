%---------------------------------------------------
% This script simulates the Bunched Phase Encoding
% which based on Fourier shift theorem
% Reference: MRM, 2006, 55:633-648
% Created by Jia Sen on 2017/04/06
%---------------------------------------------------
clc;clear;close all;

%% Simulate a finite-support 1D image
N = 256;                      % image matrix size
delta_x = 1;                  % image resolution
x = -N/2:delta_x:N/2-1;       % image pixel positions
im = zeros(N,1);
im(N/2-100:1:N/2+100-1) = 1;  
im(N/2-50:1:N/2+50-1) = 1.5;  % simulate a 1D finite-support imaged object
figure;plot(x,im);axis([min(x),max(x),1.5*min(im),1.5*max(im)]);title('original 1D finite support signal')


%% Simulate a 10x oversampled kspace
delta_k = delta_x/N;                % Nyquist sampling interveal
delta_k_os = delta_k/10;            % oversampling
k_os = 0:delta_k_os:1-delta_k_os;   % oversampled k-space grid
kspace_os = zeros(size(k_os));      % oversampled k-space

for iter_k_os = 1:1:length(k_os)    % simulate k-space sampling
    kspace_os(iter_k_os) = sum(im'.*exp(-1i*2*pi*x*k_os(iter_k_os)));
end
figure;
subplot(2,1,1);plot(k_os,abs(fftshift(kspace_os)),'k-');axis([min(k_os),max(k_os),0,1.5*max(abs(kspace_os))]);title('10x oversampled kspace: abs');
subplot(2,1,2);plot(k_os,angle(fftshift(kspace_os)),'k.');axis([min(k_os),max(k_os),-4,4]);title('10x oversampled kspace: phase');


%% Simulate Nyquist sampled kspace
delta_k = delta_x/N;                % Nyquist sampling interveal
k_nyquist = 0:delta_k:1-delta_k;    % Note, kspace range from [0,1], not centered [-0.5,0.5]
kspace_nyquist = zeros(size(k_nyquist));  % Nyquist sampled k-space
for iter_k = 1:1:length(k_nyquist)  % simulate k-space sampling via DFT of image
    kspace_nyquist(iter_k) = sum(im'.*exp(-1i*2*pi*x*k_nyquist(iter_k)));
end
figure;
subplot(2,1,1);plot(k_nyquist,abs(fftshift(kspace_nyquist)));axis([min(k_nyquist),max(k_nyquist),0,1.5*max(abs(kspace_nyquist))]);title('Nyquist sampled kspace: abs');
subplot(2,1,2);plot(k_nyquist,angle(fftshift(kspace_nyquist)),'k.');axis([min(k_nyquist),max(k_nyquist),-4,4]);title('Nyquist sampled kspace: phase');

im_recon_nyquist = zeros(size(im));
for iter_im = 1:1:length(im);       % inverse Fourier transform of kspace gives the image
    im_recon_nyquist(iter_im) = (1/length(kspace_nyquist))*sum(kspace_nyquist.*exp(1i*2*pi*k_nyquist*x(iter_im)));
end
figure;
subplot(2,1,1);plot(x,abs(im_recon_nyquist));axis([min(x),max(x),-0.5,2]);title('abs of Recon signal from Nyquist sampled data')
subplot(2,1,2);plot(x,angle(im_recon_nyquist),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of Recon signal from Nyquist sampled data')
knum_plot = 25;
figure;plot(k_nyquist(1:1:2*knum_plot),abs(kspace_nyquist(1:1:2*knum_plot)),'k*');
hold on;plot(k_os(1:1:10*knum_plot),abs(kspace_os(1:1:10*knum_plot)),'k-','linewidth',1);
legend('nyquist sampled kspace','10x oversampled kspace');
axis([k_os(1),k_os(10*knum_plot),0,1.2*max(abs(kspace_os(1:1:10*knum_plot)))]);
title('Nyquist sampled kspace')

%% simulate 2x downsampling
ds_factor = 2;
k_ds = 0:ds_factor*delta_k:1-ds_factor*delta_k;
kspace_ds = zeros(size(k_ds));
for iter_k = 1:1:length(k_ds)      % simulate k-space sampling
    kspace_ds(iter_k) = sum(im'.*exp(-1i*2*pi*x*k_ds(iter_k)));
end
figure;
subplot(2,1,1);plot(k_ds,abs(fftshift(kspace_ds)),'-k.');axis([min(k_ds),max(k_ds),0,1.5*max(abs(kspace_ds))]);title('2x downsampled kspace: abs');
subplot(2,1,2);plot(k_ds,angle(fftshift(kspace_ds)),'k.');axis([min(k_ds),max(k_ds),-3.5,3.5]);title('2x downsampled kspace: phase');

knum_plot = 25;
figure;plot(k_ds(1:1:knum_plot),abs(kspace_ds(1:1:knum_plot)),'ko');
hold on;plot(k_nyquist(1:1:2*knum_plot),abs(kspace_nyquist(1:1:2*knum_plot)),'k*');
hold on;plot(k_os(1:1:10*knum_plot),abs(kspace_os(1:1:10*knum_plot)),'k-','linewidth',1);
legend('2x downsampled kspace','nyquist sampled kspace','10x oversampled kspace');
axis([k_os(1),k_os(10*knum_plot),0,1.2*max(abs(kspace_os(1:1:10*knum_plot)))]);
title('2x downsampled kspace')

im_recon_aliase = zeros(size(im));
for iter_im = 1:1:length(im);      % inverse Fourier transform of kspace gives the aliased image
    im_recon_aliase(iter_im) = (0.5/length(kspace_ds))*sum(kspace_ds.*exp(1i*2*pi*k_ds*x(iter_im)));
end
figure;
subplot(2,1,1);plot(x,abs(im_recon_aliase));axis([min(x),max(x),-0.5,2]);title('abs of Recon signal from downsampled data')
subplot(2,1,2);plot(x,angle(im_recon_aliase),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of Recon signal from downsampled data')



%% simulate bunched 2x downsampling
alpha = 0.2;                        % distance between bunched samples is alpha*delta_k
k_shift = alpha*delta_k;            % the k-space shift should be smaller than Nyquist sampling interveral        
k_ds_shift = 0+k_shift:2*delta_k:1; % k-space sample positions
kspace_ds_shift = zeros(size(k_ds_shift));
for iter_k = 1:1:length(k_ds_shift)     % simulate k-space sampling
    kspace_ds_shift(iter_k) = sum(im'.*exp(-1i*2*pi*x*k_ds_shift(iter_k)));
end
figure;
subplot(2,1,1);plot(k_ds_shift,abs(fftshift(kspace_ds_shift)),'-k.');axis([min(k_ds_shift),max(k_ds_shift),0,1.5*max(abs(kspace_ds_shift))]);title('2x downsampled&shifted kspace: abs');
subplot(2,1,2);plot(k_ds_shift,angle(fftshift(kspace_ds_shift)),'k.');axis([min(k_ds_shift),max(k_ds_shift),-4,4]);title('2x downsampled&shifted kspace: phase');

knum_plot = 25;
figure;plot(k_ds_shift(1:1:knum_plot),abs(kspace_ds_shift(1:1:knum_plot)),'ks');
hold on;plot(k_ds(1:1:knum_plot),abs(kspace_ds(1:1:knum_plot)),'ko');
hold on;plot(k_nyquist(1:1:2*knum_plot),abs(kspace_nyquist(1:1:2*knum_plot)),'k*');
hold on;plot(k_os(1:1:10*knum_plot),abs(kspace_os(1:1:10*knum_plot)),'k-','linewidth',1);
legend('2x downsampled and shift kspace','2x downsampled kspace','nyquist sampled kspace','10x oversampled kspace');
axis([k_os(1),k_os(10*knum_plot),0,1.2*max(abs(kspace_os(1:1:10*knum_plot)))]);
title('2x bunched downsampled kspace')


im_recon_aliase_phasemod = zeros(size(im));
for iter_im = 1:1:length(im);          % inverse Fourier transform gives the aliased image with phase modulation
    im_recon_aliase_phasemod(iter_im) = (0.5/length(kspace_ds_shift))*sum(kspace_ds_shift.*exp(1i*2*pi*k_ds_shift*x(iter_im)));
end
figure;
subplot(2,1,1);plot(x,abs(im_recon_aliase_phasemod));axis([min(x),max(x),-0.5,2]);title('abs of Recon signal from downsampled ans shifted data')
subplot(2,1,2);plot(x,angle(im_recon_aliase_phasemod),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of Recon signal from downsampled and shifted data')

%% image reconstruction via BPE matrix inversion
% given im_recon_aliase and im_recon_aliase_phasemod and alpha, recon a
% im_recon_inv without aliasing artifacts via matrix inversion
im_recon_inv = zeros(size(im));
for iter_im = 1:1:N/2
    A = 0.5*[1,1;1,exp(-1i*alpha*pi)];
    b = [im_recon_aliase(iter_im);im_recon_aliase_phasemod(iter_im)];
    x_temp = A\b;
    im_recon_inv(iter_im) = x_temp(1);
end
for iter_im = N/2+1:1:N
    A = 0.5*[1,1;1,exp(1i*alpha*pi)];
    b = [im_recon_aliase(iter_im);im_recon_aliase_phasemod(iter_im)];
    x_temp = A\b;
    im_recon_inv(iter_im) = x_temp(1);
end
figure;
subplot(2,1,1);plot(x,abs(im_recon_inv));axis([min(x),max(x),-0.5,2]);title('abs of BPE Recon signal')
subplot(2,1,2);plot(x,angle(im_recon_inv),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of BPE Recon signal');

recon_error = im - abs(im_recon_inv); figure;plot(x,recon_error);title('BPE recon error')
