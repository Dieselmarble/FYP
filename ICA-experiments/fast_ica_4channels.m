% --- Bi-Orthogonal wavelet transform (Wavelab is mandatory)
% --- Denoised Sources are direct outputs of the blind-GMCA algorithm
% --- Raw Sources are obtained by applying the pseudo-inverse of the estimated mixing matrix to the data

close all; clear all; clc;

SNR_db = 30;  %-- SNR in dB

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

Xw = A*[reshape(wc1,1,numel(wc1));reshape(wc2,1,numel(wc2));reshape(wc3,1,numel(wc3));reshape(wc4,1,numel(wc4))]; %--- input data for fgmca

Xw = Xw + std2(Xw)*10^(-SNR_db/20)*randn(size(Xw));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[piA,S] = fgmca(Xw,4,200,3); %optimal parameter (800,5) for soft thresholding ; (200,3) for hard thresholding

Sw =  piA*Xw;

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

figure
subplot(221)
imnb(X1)
title('Mixture 1')
subplot(222)
imnb(X2)
title('Mixture 2')
subplot(223)
imnb(X3)
title('Mixture 3')
subplot(224)
imnb(X4)
title('Mixture 4')

% rearrange permuataion order of estimated sources
s_collect = {im1;im2;im3;im4};
esti_collect = {Sdn1;Sdn2;Sdn3;Sdn4};
s = [Sdn1(1,:);Sdn2(1,:);Sdn3(1,:);Sdn4(1,:)];
% s = [im1(1,:);im2(1,:);im3(1,:);im4(1,:)];
for j = 1:4
    cor = zeros(1,4);
    for i = 1:4
        cor(i) = abs(corr2(esti_collect{j},s_collect{i}));
    end
    [M,I] = max(cor);
    temp = esti_collect{j};
    esti_collect{j} = esti_collect{I};
    esti_collect{I} = temp;
end

Sdn1 = esti_collect{1};Sdn2 = esti_collect{2};Sdn3 = esti_collect{3};Sdn4 = esti_collect{4};
s_ = [Sdn1(1,:);Sdn2(1,:);Sdn3(1,:);Sdn4(1,:)];

% P = floor(abs(s_*pinv(s))); % get permutation matrix p
P = abs(s_*pinv(s));
% P = s_*pinv(piA*A*s);
figure
subplot(221)
imnb(Sdn1)
title('Thresholded Source - GMCA Output - 1')
subplot(222)
imnb(Sdn2)
title('Thresholded Source - GMCA Output - 2')
subplot(223)
imnb(Sdn3)
title('Thresholded Source - GMCA Output - 3')
subplot(224)
imnb(Sdn4)
title('Thresholded Source - GMCA Output - 4')

figure
subplot(221)
imnb(im1)
title('Raw Source - 1')
subplot(222)
imnb(im2)
title('Raw Source - 2')
subplot(223)
imnb(im3)
title('Raw Source - 3')
subplot(224)
imnb(im4)
title('Raw Source - 4')

% calculating psnr between estimated and original sources
psnr1 = psnr(im1,Sdn1);
psnr2 = psnr(im2,Sdn2);
psnr3 = psnr(im3,Sdn3);
psnr4 = psnr(im4,Sdn4);
% calculating correlation between estimated and original sources
R1 = abs(corr2(im1,Sdn1));
R2 = abs(corr2(im2,Sdn2));
R3 = abs(corr2(im3,Sdn3));
R4 = abs(corr2(im4,Sdn4));

cri = mmc(A,piA,P);

% [MER,perm] = bss_eval_mix(pinv(piA),A);