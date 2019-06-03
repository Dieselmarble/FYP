% --- Bi-Orthogonal wavelet transform (Wavelab is mandatory)
% --- Denoised Sources are direct outputs of the blind-GMCA algorithm
% --- Raw Sources are obtained by applying the pseudo-inverse of the estimated mixing matrix to the data

close all; clear all; clc;

SNR_db = 30;  %-- SNR in dB

iS = 256;
nTRs = 256;
nComps = 2;
nMixtures = 4;
pixles = 65536;
%sigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im1 = double(imread('boat256.png'));
im2 = double(imread('paraty256_2.tif'));
im3 = double(imread('barbara256.png'));
im4 = double(imread('pakhawaj.tif'));

im1 = im1 - mean(reshape(im1,1,numel(im1)))*ones(size(im1));
im2 = im2 - mean(reshape(im2,1,numel(im2)))*ones(size(im2));
im3 = im3 - mean(reshape(im3,1,numel(im1)))*ones(size(im1));
im4 = im4 - mean(reshape(im4,1,numel(im2)))*ones(size(im2));

% Observed image (mixtures)
Am = randn(nMixtures,nComps);
x = Am*[reshape(im1,1,pixles);reshape(im2,1,pixles)]%;reshape(im3,1,pixles);reshape(im4,1,pixles)]; 
x = reshape(x.',1,[]);
% number of PSNR coefficients to enumerate
list_PSNR = [30:5:40];
nIter = length(list_PSNR);
% figure settings
R_all = [];
noise = randn(size(x));
Y = [];
cri_all = [];
for iter = 1:nIter
    % define the PSNR peak signal to noise ratio imshow(reshape(ans1(:,2),256,256))
    PSNR  = list_PSNR(iter);
    sigma = std(x(:))*10^(-PSNR/20);
    y = x + sigma*noise; 
    % perform fastICA 
    tempy = reshape(y,[65536,nMixtures]);
    [icasig,A,W] = fastica(tempy','numOfIC',nComps, 'verbose','off');
    % get the correct permuation matrix 
    [MER,perm]=bss_eval_mix(Am,A);
    P = zeros(size(perm,1),size(perm,1));
    for i = 1:size(perm,1)
        P(i,perm(i)) = 1;
    end
    % mxing matrix criterion
    cri =  mmc_ica(Am,W,P);
    cri_all = [cri_all,cri];
    % ---------------------------- %
    % estimated sources
    icasig = P*icasig;
    for iter=1:nComps
        C{iter}=reshape(icasig(iter,:),[iS iS]);
    end;
%     imshow(C{1});
    R1 = corr2(normalize(im1),normalize(C{1}));
    R2 = corr2(normalize(im2),normalize(C{2}));
    R1 = abs(R1);
    R2 = abs(R2);
    R_all = [R_all;R1];
%     fprintf('correlation between MCA contour and original contour: %d \n', R1);
%     fprintf('correlation between MCA texture and original texture: %d \n', R2);
    
end;

%plotting
figure;
plot(list_PSNR,R_all,'LineWidth',2);
xlabel('PSNR value') 
ylabel('Correlation') 
figure;
plot(list_PSNR,cri_all,'LineWidth',2);
xlabel('PSNR value') 
ylabel('Mixing Matrix Criterion') 