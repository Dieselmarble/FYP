clear all; close all;
A = imread('lena.jpg'); 
% A=OpenBitmap('6.tif');
k=64;
l=64;
S_CI=DvdBptSubBp(A,k,l); %split into patches
% normalise
% S_CI = normalize(S_CI);
% whitening
% S_CI = whiten(S_CI);
%-----------------------------------------------------------------
[icasig, A_matrix, W_matrix] =fastica(S_CI); %independent signals; mixing matrix; sepration matrix
% -----frequency spectrum ------
% figure;
% for i = 1:49
%     Fs = 2000;
%     t = 0:1/Fs:1-1/Fs;
%     N = length(icasig(i,:));
%     xdft = fft(icasig(i,:));
%     xdft = xdft(1:N/2+1);
%     psdx = (1/(Fs*N)) * abs(xdft).^2;
%     psdx(2:end-1) = 2*psdx(2:end-1);
%     freq = 0:Fs/length(icasig(i,:)):Fs/2;
%     hold on
%     plot(freq,10*log10(psdx));
% end
% hold off
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% ------- histogram ---------------
figure;
histogram(A_matrix);
% A_matrix is 'ICA coefficients' which has super guassian distribution
% icasig is the basis funcions
% show all ica feature basis
% icasig(49,:)=0.001;
% A_matrix(:,49)=0.001;
% W_matrix(49,:)=0.001;
% figure;
% for i = 1:49
%     subplot(7,7,i);
%     imshow(reshape(icasig(i,:),64,64));
% end
%---Step1------------------------------------------------------------------
mark = imread('changsha.bmp');
b_CI=DvdBptSubBp(mark,k,l);
[ica_mark, Ab_matrix, Wb_matrix] =fastica(b_CI);
W=WMToV(mark); %wartermark.bmp

%--------------------------------------------------------------------------
index=max_cov(icasig); %find maximum covariance(energy) in decomposed signal
S_W=[icasig(index,:); W];
SS = [1 1;10 0];
%-----Step3:Y=A*S_W--------------------------------------------------------
Y=SS*double(S_W); %Y is a 2X2 matrix, linear combination of two observations
y1=Y(1,:);
y2=Y(2,:); 
%--------------------------------------------------------------------------
icasig(index,:)=y1;
S_CI=inv(W_matrix)*icasig; % Remix the signals
A_embed = SubBPtBP(A,k,l,S_CI); % Restore watermarked image
%--------------------------------------------------------------------------
figure;
imshow(A_embed);
psnr_cover = 10*log10(psnr(A,A_embed,480,480));