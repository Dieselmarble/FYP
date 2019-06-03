% --- Bi-Orthogonal wavelet transform (Wavelab is mandatory)
% --- Denoised Sources are direct outputs of the blind-GMCA algorithm
% --- Raw Sources are obtained by applying the pseudo-inverse of the estimated mixing matrix to the data

close all; clear all; clc;
nc = 10;  %-- Number of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im1 = double(imread('boat256.png'));
im2 = double(imread('paraty256_2.tif'));
im3 = double(imread('barbara256.png'));
im4 = double(imread('pakhawaj.tif'));

im1 = im1 - mean(reshape(im1,1,numel(im1)))*ones(size(im1));
im2 = im2 - mean(reshape(im2,1,numel(im2)))*ones(size(im2));
im3 = im3 - mean(reshape(im3,1,numel(im1)))*ones(size(im1));
im4 = im4 - mean(reshape(im4,1,numel(im2)))*ones(size(im2));


qmf = MakeONFilter('Battle',5);
wc1 = FWT2_PO(double(im1),1,qmf);
wc2 = FWT2_PO(double(im2),1,qmf);
wc3 = FWT2_PO(double(im3),1,qmf);
wc4 = FWT2_PO(double(im4),1,qmf);

A = randn(nc,4);

Xw_ = A*[reshape(wc1,1,numel(wc1));reshape(wc2,1,numel(wc2));reshape(wc3,1,numel(wc3));reshape(wc4,1,numel(wc4))]; %--- input data for fgmca

SNR_list = [3:5:20];
cri_all = [];
R_all = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for noise = 1:length(SNR_list)
    SNR_db = SNR_list(noise);
    Xw = Xw_ + std2(Xw_)*10^(-SNR_db/20)*randn(size(Xw_));
    % use fast gmca algorithm
    [piA,S] = fgmca(Xw,4,200,3); %optimal parameter (800,5) for soft thresholding ; (200,3) for hard thresholding
    [MER,perm] = bss_eval_mix(pinv(piA),A);
    P = zeros(size(perm,1),size(perm,1));
    for i = 1:size(perm,1)
        P(i,perm(i)) = 1;
    end 
    % P = abs(s_*pinv(s)); % get permutation matrix P
    Sw =  piA*Xw;
    Sw = P*Sw;
    % choose the firt 4 mixtures
    wwc1 = reshape(Xw(1,:),length(im1),length(im1)); 
    wwc2 = reshape(Xw(2,:),length(im1),length(im1));
    wwc3 = reshape(Xw(3,:),length(im1),length(im1));
    wwc4 = reshape(Xw(4,:),length(im1),length(im1));

    X1 = IWT2_PO(wwc1,1,qmf); % mixtures
    X2 = IWT2_PO(wwc2,1,qmf);
    X3 = IWT2_PO(wwc3,1,qmf);
    X4 = IWT2_PO(wwc4,1,qmf);

    wwc1 = reshape(S(1,:),length(im1),length(im1)); % wwc1 is estimated sources
    wwc2 = reshape(S(2,:),length(im1),length(im1));
    wwc3 = reshape(S(3,:),length(im1),length(im1));
    wwc4 = reshape(S(4,:),length(im1),length(im1));

    Sdn1 = IWT2_PO(wwc1,1,qmf); %GMCA outputs
    Sdn2 = IWT2_PO(wwc2,1,qmf);
    Sdn3 = IWT2_PO(wwc3,1,qmf);
    Sdn4 = IWT2_PO(wwc4,1,qmf);

    wwc1 = reshape(Sw(1,:),length(im1),length(im1));
    wwc2 = reshape(Sw(2,:),length(im1),length(im1));
    wwc3 = reshape(Sw(3,:),length(im1),length(im1));
    wwc4 = reshape(Sw(4,:),length(im1),length(im1));

    S1 = IWT2_PO(wwc1,1,qmf);
    S2 = IWT2_PO(wwc2,1,qmf);
    S3 = IWT2_PO(wwc3,1,qmf);
    S4 = IWT2_PO(wwc4,1,qmf);
    
    % calculating correlation between estimated and original sources
    R1 = abs(corr2(normalize(im1),normalize(S1)));
    R2 = abs(corr2(normalize(im2),normalize(S2)));
    R3 = abs(corr2(normalize(im3),normalize(S3)));
    R4 = abs(corr2(normalize(im4),normalize(S4)));
    R_all = [R_all;R1];
    cri = mmc(A,piA,P); 
    cri_all = [cri_all;cri];
end
% --------- plot obervations ---------------%
figure;
plot(SNR_list,R_all);
xlabel('PSNR value') 
ylabel('Correlation Coefficient') 
figure;
plot(SNR_list,cri_all);
xlabel('PSNR value') 
ylabel('Mixing Matrix Criterion') 