clear all; close all;
A = imread('lena.jpg');
% A=OpenBitmap('6.tif');
k=64;
l=64; %32
S_CI=DvdBptSubBp(A,k,l); %split into patches
[icasig, A_matrix, W_matrix] =fastica(S_CI); %independent signals; mixing matrix; sepration matrix
%---Step1------------------------------------------------------------------
mark = imread('changsha.bmp');
W=WMToV(mark); %wartermark.bmp
%--------------------------------------------------------------------------
index=max_cov(icasig); %find maximum covariance(energy) in decomposed signal
S_W=[icasig(index,:); W];
SS = [0.5 0.5; 0.1 0.9];
%--------------------------------------------------------------------------

%-----Step3:Y=A*S_W--------------------------------------------------------
Y=SS*S_W; %Y is a 2X2 matrix, linear combination of two observations
y1=Y(1,:);
y2=Y(2,:); 
%--------------------------------------------------------------------------
icasig(index,:)=y1; %embedded watermark in max covariance
S_CI=inv(W_matrix)*icasig; % Remix the signals
SubBPtBP(A,k,l,S_CI); % Restore watermarked image
%--------------------------------------------------------------------------
figure(1);
imshow(A);






