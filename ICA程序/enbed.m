clear all;
A=OpenBitmap('6.tif');
k=25;
l=32;
S_CI=DvdBptSubBp(A,k,l); %split into patches
[icasig, A_matrix, W_matrix] =fastica(S_CI);
%---Step1------------------------------------------------------------------
mark = imread('wartermark.bmp');
W=WMToV(mark); %wartermark.bmp
%--------------------------------------------------------------------------
index=max_cov(icasig); %find maximum covariance(energy) in decomposed signal
S_W=[icasig(index,:); W];
SS = [0.8 0.1; 0.9 0.2];
%--------------------------------------------------------------------------

%-----Step3:Y=A*S_W--------------------------------------------------------
Y=SS*S_W;
y1=Y(1,:);
y2=Y(2,:);
%--------------------------------------------------------------------------
icasig(index,:)=y1;
S_CI=inv(W_matrix)*icasig;%
SubBPtBP(A,k,l,S_CI);
%--------------------------------------------------------------------------
figure(1);
imshow(A);






