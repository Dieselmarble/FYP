%%%%%%%%%%%%%%% Boy image %%%%%%%%%%%%%%%%%%%%%%%%
clear all; clear all; clc;
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
img1 = vecOrig1 + vecOrig2;
img2 = vecOrig1 + vecOrig2;
x = [img1,img2];
% number of PSNR coefficients to enumerate
nIter = 7;
list_PSNR = [0:5:30];
% figure settings
R_all = [];
noise = randn(size(x));
Y = [];
for iter = 1:nIter,
    % define the PSNR peak signal to noise ratio
    PSNR  = list_PSNR(iter);
    sigma = std(x(:))*10^(-PSNR/20);
    y = x + sigma*noise; 
    Y = [Y;y];
    % perform fastICA 
    [icasig,A,W] = fastica (reshape(y,[65536,nMixtures]).', 'numOfIC', nComps, 'verbose','off');
    for iter=1:nComps,
        C{iter}=reshape(icasig(iter,:),[iS iS]);
    end;
    R1 = corr2(imgc,C{1});
    R2 = corr2(imgt,C{2});
    R1 = abs(R1);
    R2 = abs(R2);
    R_all = [R_all;R1];
    fprintf('correlation between MCA contour and original contour: %d \n', R1);
    fprintf('correlation between MCA texture and original texture: %d \n', R2);
end;

%plotting
plot(list_PSNR,R_all);

