%---------------------------------------------------
% This script simulates the Bunched Phase Encoding
% which based on Fourier shift theorem
% Generalized version. 
% ----------------------
% Algorithm parameters:
% ds_factor: downsampling factor
% m: bunched sample size
% alpha: <1, distance between adjacent bunched samples
% ----------------------
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
im(N/2-30:1:N/2-1) = 1.5; 
im(N/2:1:N/2+15) = 3;          % simulate a 1D finite-support imaged object

% im = sin(2*pi*(x+N/2)/(2*N));
% im = reshape(im,[N,1]);       % simulate a 1D finite-support imaged object

figure;plot(x,im);axis([min(x),max(x),-0.5,1.5*max(im)]);title('original 1D finite support signal')


%% Simulate an oversampled kspace
delta_k = delta_x/N;                % Nyquist sampling interveal
os_factor = 16;
delta_k_os = delta_k/os_factor;            % oversampling
k_os = 0:delta_k_os:1-delta_k_os;   % oversampled k-space grid
kspace_os = zeros(size(k_os));      % oversampled k-space
for iter_k_os = 1:1:length(k_os)    % simulate k-space sampling
    kspace_os(iter_k_os) = sum(im'.*exp(-1i*2*pi*x*k_os(iter_k_os)));
end

figure;
subplot(2,1,1);plot(k_os,abs(fftshift(kspace_os)),'k-');axis([min(k_os),max(k_os),0,1.5*max(abs(kspace_os))]);title('oversampled kspace: abs');
subplot(2,1,2);plot(k_os,angle(fftshift(kspace_os)),'k.');axis([min(k_os),max(k_os),-4,4]);title('oversampled kspace: phase');


%% Simulate Nyquist sampled kspace
delta_k = delta_x/N;                % Nyquist sampling interveal
k_nyquist = 0:delta_k:1-delta_k;    % Note, kspace range from [0,1], not centered [-0.5,0.5]
kspace_nyquist = zeros(size(k_nyquist));  % Nyquist sampled k-space
for iter_k = 1:1:length(k_nyquist)  % simulate k-space sampling via DFT of image
    kspace_nyquist(iter_k) = sum(im'.*exp(-1i*2*pi*x*k_nyquist(iter_k)));
end

im_recon_nyquist = zeros(size(im));
for iter_im = 1:1:length(im);       % inverse Fourier transform of kspace gives the image
    im_recon_nyquist(iter_im) = (1/N)*sum(kspace_nyquist.*exp(1i*2*pi*k_nyquist*x(iter_im)));
end

figure;
subplot(2,1,1);plot(k_nyquist,abs(fftshift(kspace_nyquist)));axis([min(k_nyquist),max(k_nyquist),0,1.5*max(abs(kspace_nyquist))]);title('Nyquist sampled kspace: abs');
subplot(2,1,2);plot(k_nyquist,angle(fftshift(kspace_nyquist)),'k.');axis([min(k_nyquist),max(k_nyquist),-4,4]);title('Nyquist sampled kspace: phase');
figure;
subplot(2,1,1);plot(x,abs(im_recon_nyquist));axis([min(x),max(x),-0.5,2]);title('abs of Recon signal from Nyquist sampled data')
subplot(2,1,2);plot(x,angle(im_recon_nyquist),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of Recon signal from Nyquist sampled data')
knum_plot = 25;
figure;plot(k_nyquist(1:1:2*knum_plot),abs(kspace_nyquist(1:1:2*knum_plot)),'ko');
hold on;plot(k_os(1:1:10*knum_plot),abs(kspace_os(1:1:10*knum_plot)),'k-','linewidth',1);
legend('nyquist sampled kspace','oversampled kspace');
axis([k_os(1),k_os(10*knum_plot),0,1.2*max(abs(kspace_os(1:1:10*knum_plot)))]);
title('Nyquist sampled kspace')


%% simulate downsampling
ds_factor = 4;     % p = ds_factor, downsampling factor
k_ds = 0:ds_factor*delta_k:1-ds_factor*delta_k;
kspace_ds = zeros(size(k_ds));
for iter_k = 1:1:length(k_ds)      % simulate k-space sampling
    kspace_ds(iter_k) = sum(im'.*exp(-1i*2*pi*x*k_ds(iter_k)));
end

im_recon_aliase = zeros(size(im));
for iter_im = 1:1:length(im);      % inverse Fourier transform of kspace gives the aliased image
    im_recon_aliase(iter_im) = (1/N)*sum(kspace_ds.*exp(1i*2*pi*k_ds*x(iter_im)));
end

figure;
subplot(2,1,1);plot(k_ds,abs(fftshift(kspace_ds)),'-k.');axis([min(k_ds),max(k_ds),0,1.5*max(abs(kspace_ds))]);title('downsampled kspace: abs');
subplot(2,1,2);plot(k_ds,angle(fftshift(kspace_ds)),'k.');axis([min(k_ds),max(k_ds),-3.5,3.5]);title('downsampled kspace: phase');


figure;
knum_plot = 16;
plot(k_ds(1:1:knum_plot),abs(kspace_ds(1:1:knum_plot)),'k*');
hold on;plot(k_nyquist(1:1:ds_factor*knum_plot),abs(kspace_nyquist(1:1:ds_factor*knum_plot)),'ko');
hold on;plot(k_os(1:1:os_factor*knum_plot),abs(kspace_os(1:1:os_factor*knum_plot)),'k-','linewidth',1);
legend('downsampled kspace','nyquist sampled kspace','oversampled kspace');
axis([k_os(1),k_os(os_factor*knum_plot),0,1.2*max(abs(kspace_os(1:1:os_factor*knum_plot)))]);
title('downsampled kspace')

figure;
subplot(2,1,1);plot(x,abs(im_recon_aliase));axis([min(x),max(x),-0.5,2]);title('abs of Recon signal from downsampled data')
subplot(2,1,2);plot(x,angle(im_recon_aliase),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of Recon signal from downsampled data')



%% simulate bunched downsampling
m = 5;           % Bunched sample size, should be larger than downsampling factor
alpha = 0.05;    % distance between adjacent bunched samples is alpha*delta_k
k_ds_bpe = zeros(N/ds_factor,m);
kspace_ds_bpe = zeros(N/ds_factor,m);
for i_alpha = 1:1:m
    k_shift = alpha*delta_k*i_alpha;               % the k-space shift should be smaller than Nyquist sampling interveral
    k_ds_shift_m = 0+k_shift:ds_factor*delta_k:1;  % k-space sample positions
    k_ds_bpe(:,i_alpha) = k_ds_shift_m;
    kspace_ds_shift_m = zeros(size(k_ds_shift_m));
    for iter_k = 1:1:length(k_ds_bpe)     % simulate k-space sampling
        kspace_ds_shift_m(iter_k) = sum(im'.*exp(-1i*2*pi*x*k_ds_shift_m(iter_k)));
    end
    kspace_ds_bpe(:,i_alpha) = kspace_ds_shift_m;
end

figure;plot(k_ds_bpe(1:1:knum_plot,2),abs(kspace_ds_bpe(1:1:knum_plot,2)),'ks');
if m >= 4 && m < 6
    hold on;plot(k_ds_bpe(1:1:knum_plot,4),abs(kspace_ds_bpe(1:1:knum_plot,4)),'ko');
    legend('BPE, m=2','BPE, m=4','BPE, m=0','oversampled kspace');
end
if m >= 6
    hold on;plot(k_ds_bpe(1:1:knum_plot,6),abs(kspace_ds_bpe(1:1:knum_plot,6)),'k+');
    legend('BPE, m=2','BPE, m=4','BPE, m=6','BPE, m=0','oversampled kspace');
end
hold on;plot(k_ds(1:1:knum_plot),abs(kspace_ds(1:1:knum_plot)),'k*');
hold on;plot(k_os(1:1:10*knum_plot),abs(kspace_os(1:1:10*knum_plot)),'k-','linewidth',1);
axis([k_os(1),k_os(10*knum_plot),0,1.2*max(abs(kspace_os(1:1:10*knum_plot)))]);
title('bunched downsampled kspace')

im_bpe_phasemod = zeros(numel(im),m);
for i_alpha = 1:1:m
    im_recon_aliase_phasemod_m = zeros(size(im));
    for iter_im = 1:1:length(im);          % inverse Fourier transform gives the aliased image with phase modulation
        im_recon_aliase_phasemod_m(iter_im) = (1/N)*sum(kspace_ds_bpe(:,i_alpha).*exp(1i*2*pi*k_ds_bpe(:,i_alpha)*x(iter_im)));
    end
    im_bpe_phasemod(:,i_alpha) = im_recon_aliase_phasemod_m;
    figure;
    subplot(2,1,1);plot(x,abs(im_bpe_phasemod(:,i_alpha)));axis([min(x),max(x),-0.5,2]);title('abs of Recon signal from downsampled and shifted data')
    subplot(2,1,2);plot(x,angle(im_bpe_phasemod(:,i_alpha)),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of Recon signal from downsampled and shifted data')
end


%% image reconstruction via BPE matrix inversion
% given im_recon_aliase and im_bpe_phasemod and alpha, recon a
% im_recon_inv without aliasing artifacts via matrix inversion
im_recon_inv = zeros(size(im));
A = zeros(m+1,ds_factor);
for i_alpha = 1:1:m+1
    for i_ds = 1:1:ds_factor
        A(i_alpha,i_ds) = 1/ds_factor * exp(-1i*2*pi*alpha*(i_alpha-1)*(i_ds - 1)/ds_factor);
    end
end
for iter_im = 1:1:N/ds_factor
    b = [im_recon_aliase(iter_im);reshape(im_bpe_phasemod(iter_im,:),[m,1])];
    x_temp = A\b;
    for i_ds = 1:1:ds_factor
        im_recon_inv(iter_im + (i_ds-1)*N/ds_factor) = x_temp(i_ds);
    end
end
im_recon_inv = im_recon_inv(1:1:numel(im));

figure;
subplot(2,1,1);plot(x,abs(im_recon_inv)); axis([min(x),max(x),-0.5,4]);title('abs of BPE Recon signal')
subplot(2,1,2);plot(x,angle(im_recon_inv),'k.');axis([min(x),max(x),-3.5,3.5]);title('phase of BPE Recon signal');

recon_error = im - abs(im_recon_inv); figure;plot(x,recon_error);title('BPE recon error')
