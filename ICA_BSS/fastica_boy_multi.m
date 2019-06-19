%%%%%%%%%%%%%%% Boy image %%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
iS = 256;
nTRs = 256;
nComps = 2;
nMixtures = 2;
%sigma;

meanValue = 0;
var_gauss = 0.3;
imgt = 0.4*double(imread('texture4.tif'));
imgc = 0.6*double(imread('boy.tif'));
vecOrig1 = imgt(:)';
vecOrig2 = imgc(:)';

% normalise original image to unity variance
vecOrig1 = normalize(vecOrig1);
vecOrig2 = normalize(vecOrig2);

% vecOrig1=vecOrig1/std(vecOrig1(:));
% vecOrig2=vecOrig2/std(vecOrig2(:));

% add Gaussian noise to image as a row vector, use this or the PSNR method
% vecOrig1 = vecOrig1 + sqrt(var_gauss)*randn(size(vecOrig1)) + meanValue;
% vecOrig2 = vecOrig2 + sqrt(var_gauss)*randn(size(vecOrig2)) + meanValue;
% y = x + sqrt(var_gauss)*randn(size(x)) + meanValue;


% Observed image (mixtures)
% img1 = vecOrig1 + vecOrig2;
% img2 = vecOrig1 + vecOrig2;
% x = [img1,img2];
% -------------------------- %
Am = randn(nMixtures,nComps);
% Am = [0.2222,  1.3325;
%     0.3223,  0.1454];
x = Am*[vecOrig1;vecOrig2];
img1 = x(1,:);
img2 = x(2,:);
x = reshape(x.',1,[]);
% number of PSNR coefficients to enumerate
list_PSNR = [3:5:40];
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
    % Y = [Y;y];
    % perform fastICA 
    tempy = reshape(y,[65536,nMixtures]);
    [icasig,A,W] = fastica (tempy.', 'numOfIC', nComps, 'verbose','off');
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
    R1 = corr2(normalize(imgt),normalize(C{1}));
    R2 = corr2(normalize(imgc),normalize(C{2}));
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