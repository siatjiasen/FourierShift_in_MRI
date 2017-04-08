clc;clear;close all
%% simulate an image with image resolution = 1
N = 128;
delta_im = 1;
FOV_im = N*delta_im;
x2D = repmat(0:1:N-1,[N,1]);    % col -> readout
y2D = repmat((0:1:N-1)',[1,N]); % row -> phase encoding
im = phantom(N);
figure;imshow(im,[]);

%% simulate the kspace with Nyquist sampling intervel delta_k = 1/FOV = 1/N
kspace = zeros(size(im));
delta_k = 1/N;
kx = (0:1:N-1)*delta_k;     % physical coordinates of kspace samples
kx2D = repmat(kx,[N,1]);    % col -> readout
ky2D = repmat(kx',[1,N]);   % row -> phase encoding
ky= kx;
for irow = 1:1:N
    for icol = 1:1:N
        kspace(irow,icol) = sum(sum(im.*exp(-1i*2*pi*(kx(icol)*x2D + y2D*ky(irow)))));
    end
end


% image reconstruction via inverse DFT from k-space samples
im_recon = zeros(size(im));
for irow = 1:1:N
    for icol = 1:1:N
        im_recon(irow,icol) = (1/N^2)*sum(sum(kspace.*exp(1i*2*pi*(kx2D*x2D(irow,icol) + y2D(irow,icol)*ky2D))));
    end
end

figure;imshow(20*log(1+abs(fftshift(kspace))),[]);title('Nyquist sampled k-space: abs');


%% simulate kspace downsampling
ds_factor = 2;
N_B = N/ds_factor;
kspace_ds = zeros(N_B,N);
ky_ds = (0:ds_factor:N-1)*delta_k;     % physical coordinates of downsampled kspace samples
kx2D_ds = repmat(kx,[N_B,1]);          % col -> readout
ky2D_ds = repmat(ky_ds',[1,N]);        % row -> phase encoding
ky= kx;
for irow = 1:1:N_B
    for icol = 1:1:N
        kspace_ds(irow,icol) = sum(sum(im.*exp(-1i*2*pi*(kx(icol)*x2D + y2D*ky_ds(irow)))));
    end
end

figure;imshow(10*log(1+abs(fftshift(kspace_ds))),[0,100]);title('Downsampled k-space: abs');


% aliased image via inverse DFT from downsampled k-space samples
im_aliase_recon = zeros(N,N);
for irow = 1:1:N
    for icol = 1:1:N
        im_aliase_recon(irow,icol) = (1/N^2)*sum(sum(kspace_ds.*exp(1i*2*pi*(kx2D_ds*x2D(irow,icol) + y2D(irow,icol)*ky2D_ds))));
    end
end

figure;imshow(abs(im_aliase_recon),[]);title('Aliasing');
figure;imshow(angle(im_aliase_recon),[-3.5,3.5]);colormap('jet');title('Phase = 0');

%% generate BPE 2D sampling pattern 
%- downsampling along PE, oversampling along RO
r = 3;    % the width of a zigzag kspace trajectory, in unit of delta_ky/os_factor
m = 4;

os_factor = 16;             % kspace oversampling factor
PE_shift = [0, 2, 3, 1];    % shift bunched samples along PE/row, 0 shift -> reference point
RO_shift = [0, 4, 8, 12];   % shift bunched samples along RO/col, 0 shift-> reference point

bpe_sampling_pattern = zeros(N*os_factor,N*os_factor);
bpe_sampling_pattern(1:ds_factor*os_factor:size(bpe_sampling_pattern),1:os_factor:size(bpe_sampling_pattern,2)) = 1;
for ishift = 1:1:length(RO_shift)
    bpe_sampling_pattern((PE_shift(ishift)+1):ds_factor*os_factor:N*os_factor,(RO_shift(ishift)+1):os_factor:N*os_factor) = 1;
end
figure;imshow(bpe_sampling_pattern(1:1:N,1:1:N));title('ZigZag Sampling');
xlabel('frequency encoding');ylabel('phase encoding')

%% simulate BPE 2D sampling
kx_shift = delta_k*RO_shift/os_factor;
ky_shift = delta_k*PE_shift/os_factor;   % physical shift

bpe_ro_num = N*m;    % bunched sample size along readout
kx_bpe = zeros(1, bpe_ro_num);
for iter_m = 1:1:m
    kx_bpe(1,iter_m:m:end) = kx + kx_shift(iter_m);
end
bpe_pe_num = N*m/ds_factor;
ky_bpe = zeros(bpe_pe_num,1);
for iter_m = 1:1:m
    ky_bpe(iter_m:m:end,1) = ky_ds + ky_shift(iter_m);
end

tic
kspace_bpe_full = zeros(N_B,N*m);
for iPE = 1:1:N_B
    for iRO = 1:1:N
        for ishift = 1:1:m
            ky_temp = ky_bpe( (iPE-1)*m + ishift);
            kx_temp = kx_bpe( (iRO-1)*m + ishift);
            kspace_bpe_full(iPE,(iRO-1)*m + ishift) = sum(sum(im.*exp(-1i*2*pi*(kx_temp*x2D + y2D*ky_temp))));
            
        end
    end
end
toc

%% BPE lead to phase modulated aliased images
im_aliase_bpe = zeros([size(im),m]);
for iter_m = 1:1:m
    kspace_bpe_m = kspace_bpe_full(:, iter_m:m:end);
    kx1D_bpe_m = kx_bpe(iter_m:m:end);
    ky1D_bpe_m = ky_bpe(iter_m:m:end);
    kx2D_bpe_m = repmat(kx1D_bpe_m,[N_B,1]);          % col -> readout
    ky2D_bpe_m = repmat(ky1D_bpe_m,[1,N]);        % row -> phase encoding
    im_aliase_recon = zeros(N,N);
    for irow = 1:1:N
        for icol = 1:1:N
            im_aliase_recon(irow,icol) = (1/N^2)*sum(sum(kspace_bpe_m.*exp(1i*2*pi*(kx2D_bpe_m*x2D(irow,icol) + y2D(irow,icol)*ky2D_bpe_m))));
        end
    end
    im_aliase_bpe(:,:,iter_m) = im_aliase_recon;
    figure;imshow(abs(im_aliase_bpe(:,:,iter_m)));title(sprintf('BPE aliasing: abs, shift = %d',iter_m-1));
    figure;imshow(angle(im_aliase_bpe(:,:,iter_m)));colormap('jet');title(sprintf('BPE aliasing: phase, shift = %d',iter_m -1));
end

%% BPE reconstruction via inversing the BPE matrix
im_recon_bpe = zeros(size(im));
for iter_PE = 1:1:N/ds_factor
    for iter_RO = 1:1:N
        A = zeros(m,ds_factor);
        for iArow = 1:1:m
            for iAcol = 1:1:ds_factor
                A(iArow,iAcol) = (1/ds_factor)*exp(-1i*2*pi*((ky_shift(iArow)/delta_k)*(iAcol - 1)/ds_factor)); % Readout shift has no effect
            end
        end
        b = reshape(squeeze(im_aliase_bpe(iter_PE,iter_RO,:)),[m,1]);
        x = A \ b;
        for ids_factor = 1:1:ds_factor
            im_recon_bpe(iter_PE + (ids_factor-1)*N/ds_factor,iter_RO) = x(ids_factor);
        end
    end
end

figure;imshow(abs(im_recon_bpe));colormap('gray');title('BPE recon: abs')
figure;imshow(angle(im_recon_bpe));colormap('jet');title('BPE recon: phase')


