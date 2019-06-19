close all ; clear all; clc;
% ---------------------------- %
cri_all_gmca = load('cri_gmca.mat');
cri_all_ksvd = load('cri_ksvd.mat');
cri_all_bksvd = load('cri_bksvd.mat');

% ---------------------------- %
psnr_all_gmca = load('psnr_gmca.mat');
psnr_all_ksvd = load('psnr_ksvd.mat');
psnr_all_bksvd = load('psnr_bksvd.mat');


% ---------------------------- %
R_all_gmca = load('R_gmca.mat');
R_all_ksvd = load('R_ksvd.mat');
R_all_bksvd = load('R_bksvd.mat');

% ---------------------------- %
list_PSNR = [3:5:40];
figure;
hold on
plot(list_PSNR,R_all_gmca.R_all,'-kd','LineWidth',2); % Solid line (default)
plot(list_PSNR,R_all_ksvd.R_all,'--ko','LineWidth',2); % Dashed line
plot(list_PSNR,R_all_bksvd.R_all,':k*','LineWidth',2); % Dotted line
hold off
xlabel('PSNR in dB') 
ylabel('Correlation') 
legend('GMCA','K-SVD','BK-SVD','Location','southeast')
% --------------------------- %
figure;
hold on;
plot(list_PSNR,cri_all_gmca.cri_all,'-kd','LineWidth',2);
plot(list_PSNR,cri_all_ksvd.cri_all,'--ko','LineWidth',2);
plot(list_PSNR,cri_all_bksvd.cri_all,':k*','LineWidth',2);
hold off;
xlabel('PSNR in dB') 
ylabel('Mixing Matrix Criterion')
legend('GMCA','K-SVD','BK-SVD','Location','northeast')
% --------------------------- %
figure;
hold on;
plot(list_PSNR,psnr_all_gmca.psnr_all,'-kd','LineWidth',2);
plot(list_PSNR,psnr_all_ksvd.psnr_all,'--ko','LineWidth',2);
plot(list_PSNR,psnr_all_bksvd.psnr_all,':k*','LineWidth',2);
hold off;
xlabel('PSNR in dB') 
ylabel('Sinal to Noise Ratio Recoved')
legend('GMCA','K-SVD','BK-SVD','Location','southeast')
% --------------------------- %
