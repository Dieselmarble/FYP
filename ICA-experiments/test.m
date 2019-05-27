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
% img1 = vecOrig1 + vecOrig2;
% img2 = vecOrig1 + vecOrig2;
% x = [img1,img2];
% -------------------------- %
Am = randn(2,2);
x = Am*[vecOrig1;vecOrig2];
img1 = x(1,:);
img2 = x(2,:);
x = reshape(x.',1,[]);
% number of PSNR coefficients to enumerate
list_PSNR = [70:5:80];
nIter = length(list_PSNR);
% figure settings
R_all = [];
noise = randn(size(x));    
% define the PSNR peak signal to noise ratio imshow(reshape(ans1(:,2),256,256))
iter = 3;
PSNR  = list_PSNR(iter);
sigma = std(x(:))*10^(-PSNR/20);
y = x + sigma*noise; 
%     Y = [Y;y];
% perform fastICA 
tempy = reshape(y,[65536,nMixtures]);
[icasig,A,W] = fastica (tempy.', 'numOfIC', nComps, 'verbose','off');