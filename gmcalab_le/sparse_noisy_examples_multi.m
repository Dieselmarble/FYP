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
% 
% im1 = im1 - mean(reshape(im1,1,numel(im1)))*ones(size(im1));
% im2 = im2 - mean(reshape(im2,1,numel(im2)))*ones(size(im2));
% im3 = im3 - mean(reshape(im3,1,numel(im1)))*ones(size(im1));
% im4 = im4 - mean(reshape(im4,1,numel(im2)))*ones(size(im2));


qmf = MakeONFilter('Coiflet',5);
wc1 = FWT2_PO(double(im1),1,qmf);
wc2 = FWT2_PO(double(im2),1,qmf);
wc3 = FWT2_PO(double(im3),1,qmf);
wc4 = FWT2_PO(double(im4),1,qmf);

A = randn(nc,4);

Xw_ = A*[reshape(wc1,1,numel(wc1));reshape(wc2,1,numel(wc2));reshape(wc3,1,numel(wc3));reshape(wc4,1,numel(wc4))]; %--- input data for fgmca

SNR_list = [3:5:40];
cri_all = [];
R_all = [];
% SDR_all = [];
% SIR_all = [];
% SAR_all = [];
psnr_all = []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for noise = 1:length(SNR_list)
    SNR_db = SNR_list(noise);
    Xw = Xw_ + std2(Xw_)*10^(-SNR_db/20)*randn(size(Xw_));
    % use fast gmca algorithm
    [piA,S,his1] = fgmca(Xw,4,100,3); %optimal parameter (800,5) for soft thresholding ; (200,3) for hard thresholding
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
    im_source1 = mat2gray(im1)*255;
    im_source_noisy1 = mat2gray(reshape(S1,size(im1)))*255;
    PSNR_source1 = 20*log10(255 * sqrt(numel(im_source1)) / norm(im_source1(:)-im_source_noisy1(:)));

    im_source2 = mat2gray(im2)*255;
    im_source_noisy2 = mat2gray(reshape(S2,size(im2)))*255;
    PSNR_source2 = 20*log10(255 * sqrt(numel(im_source2)) / norm(im_source2(:)-im_source_noisy2(:)));

    im_source3 = mat2gray(im3)*255;
    im_source_noisy3 = mat2gray(reshape(S3,size(im3)))*255;
    PSNR_source3 = 20*log10(255 * sqrt(numel(im_source3)) / norm(im_source3(:)-im_source_noisy3(:)));

    im_source4 = mat2gray(im4)*255;
    im_source_noisy4 = mat2gray(reshape(S4,size(im4)))*255;
    PSNR_source4 = 20*log10(255 * sqrt(numel(im_source4)) / norm(im_source4(:)-im_source_noisy4(:)));

    avg_PNSR_source = (PSNR_source1 + PSNR_source2 + PSNR_source3 + PSNR_source4)/4;

    R1 = abs(corr2(im1,im_source_noisy1));
    R2 = abs(corr2(im2,im_source_noisy2));
    R3 = abs(corr2(im3,im_source_noisy3));
    R4 = abs(corr2(im4,im_source_noisy4));
    
    avg_R = (R1+R2+R3+R4)/4;
    R_all = [R_all;avg_R];
    psnr_all = [psnr_all;avg_PNSR_source];
    cri = mmc(A,piA,P); 
    cri_all = [cri_all;cri];
    
%     estimated = [Sdn1(:)';Sdn2(:)';Sdn3(:)';Sdn4(:)'];
%     true = [im1(:)';im2(:)';im3(:)';im4(:)'];
%     [SDR,SIR,SAR,perm] = bss_eval_sources(estimated,true);
%     SDR_all = [SDR_all;SDR];
%     SIR_all = [SIR_all;SIR];
%     SAR_all = [SAR_all;SAR];

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