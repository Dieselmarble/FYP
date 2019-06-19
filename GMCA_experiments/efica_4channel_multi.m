% --- Bi-Orthogonal wavelet transform (Wavelab is mandatory)
% --- Denoised Sources are direct outputs of the blind-GMCA algorithm
% --- Raw Sources are obtained by applying the pseudo-inverse of the estimated mixing matrix to the data

close all; clear all; clc;

SNR_db = 30;  %-- SNR in dB

iS = 256;
nTRs = 256;
nComps = 4;
nMixtures = 10;
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
x = Am*[reshape(im1,1,pixles);reshape(im2,1,pixles);reshape(im3,1,pixles);reshape(im4,1,pixles)]; 
x = reshape(x.',1,[]);
% number of PSNR coefficients to enumerate
list_PSNR = [3:5:20];
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
%     [icasig,A,W] = fastica(tempy','numOfIC',nComps, 'verbose','off');
%     W = jader(tempy',nComps);
    W = efica(tempy',nComps);
    % get the correct permuation matrix 
    A = pinv(W);
    [vars,keys] = sort(abs(sum(A)),'descend');
    A = A(:,keys(1:4));
%     W = W(keys(1:4),:);
    icasig = pinv(A)*tempy';
%     icasig = icasig(keys(1:4),:);
    % ------------------------ %
    s_collect = {im1;im2;im3;im4};
    for iter=1:nComps
        C{iter}=reshape(icasig(iter,:),[iS iS]);
    end;
    esti_collect = {C{1};C{2};C{3};C{4}};
    s = [C{1}(1,:);C{2}(1,:);C{3}(1,:);C{4}(1,:)];
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
    P = abs(s_*pinv(s));
    % ------------------------ %
    [MER,perm]=bss_eval_mix(Am,A);
    P1 = zeros(size(perm,1),size(perm,1));
    for i = 1:size(perm,1)
        P1(i,perm(i)) = 1;
    end
    % mxing matrix criterion
    cri =  mmc(Am,pinv(A),P);
    cri_all = [cri_all,cri];
    % ---------------------------- %
    % estimated sources
    icasig = P*icasig;
    for iter=1:nComps
        C{iter}=reshape(icasig(iter,:),[iS iS]);
    end;
%     imshow(C{1});
    im_source1 = mat2gray(im1)*255;
    im_source2 = mat2gray(im2)*255;
    im_source3 = mat2gray(im3)*255;
    im_source4 = mat2gray(im4)*255;
    im_source_noisy1 =  mat2gray(reshape(C{1},size(im1)))*255;
    im_source_noisy2 =  mat2gray(reshape(C{2},size(im2)))*255;
    im_source_noisy3 =  mat2gray(reshape(C{3},size(im3)))*255;
    im_source_noisy4 =  mat2gray(reshape(C{4},size(im4)))*255;
    R1 = abs(corr2(im1, im_source_noisy1));
    R2 = abs(corr2(im2, im_source_noisy2));
    R3 = abs(corr2(im3, im_source_noisy3));
    R4 = abs(corr2(im4, im_source_noisy4));
    avg_R = (R1+R2+R3+R4)/4;
    R_all = [R_all;avg_R];
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