%%%%%%%%%%%%%%% Boy image %%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;
iS = 256;
nTRs = 256;
nComps = 2;
%sigma;

meanValue = 0;
var_gauss = 0.3;
imgt = double(imread('texture4.tif'));
imgc = double(imread('boy.tif'));
vecOrig1 = imgt(:)';
vecOrig2 = imgc(:)';

% add Gaussian noise to image as a row vector, use this or the PSNR method
vecOrig1 = normalize(vecOrig1);
vecOrig2 = normalize(vecOrig2);
% vecOrig1 = vecOrig1 + sqrt(var_gauss)*randn(size(vecOrig1)) + meanValue;
% vecOrig2 = vecOrig2 + sqrt(var_gauss)*randn(size(vecOrig2)) + meanValue;

% Observed image (mixtures)
img1 = 2*vecOrig1 + 3*vecOrig2;
img2 = 1*vecOrig1 + 1*vecOrig2;
x = [img1,img2];
% define the PSNR peak signal to noise ratio
PSNR  = 10;
sigma = std(x(:))*10^(-PSNR/20);
noise = randn(size(x));
y = x;% + sigma*noise; 

% used by variance method, not PSNR
% y = x + sqrt(var_gauss)*randn(size(x)) + meanValue;

[icasig,A,W] = fastica (reshape(y,[65536,2]).', 'numOfIC', nComps);
%               Reduce dimension to 10, and estimate only 3
%               independent components.

for iter=1:nComps,
    C{iter}=reshape(icasig(iter,:),[iS iS]);
end;

%Display image with scaled color
%image(W,'CDataMapping','scaled');
figure;
imagesc(imgc);
colormap('gray');
axis off;
figure;
imagesc(imgt);
colormap('gray');
axis off;
figure;
imagesc(reshape(img1,256,256));
colormap('gray');
axis off;
figure;
imagesc(reshape(img2,256,256));
colormap('gray');
axis off;
figure;
imagesc(C{1});
colormap('gray');
axis off;
figure;
imagesc(C{2});
colormap('gray');
axis off;
%display the final images in gray scale
%mean square error 
err = immse(C{1},imgc);
R1 = corr2(imgc,C{1});
fprintf('correlation between MCA contour and original contour: %d \n' , R1);
R2 = corr2(imgt,C{2});
fprintf('correlation between MCA texture and original texture: %d \n', R2);

