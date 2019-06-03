close all ; clear all; clc;
% ---------------------------- %
cri_all_ica = load('cri_all_ica3.mat');
cri_all_efica = load('cri_all_efica.mat');
cri_all_jader = load('cri_all_jader3.mat');
cri_all_mmca = load('cri_all_mmca.mat');

% ---------------------------- %
R_all_ica = load('R_all_ica.mat');
R_all_efica = load('R_all_efica.mat');
R_all_jader = load('R_all_jader.mat');
R_all_mca = load('R_all_MCA.mat');
R_all_mmca = load('R_all_mmca1.mat');
% ---------------------------- %
list_PSNR = [3:5:40];
figure;
hold on
plot(list_PSNR,R_all_ica.R_all,'-kd','LineWidth',2); % Solid line (default)
plot(list_PSNR,R_all_efica.R_all,'--ko','LineWidth',2); % Dashed line
plot(list_PSNR,R_all_jader.R_all,':k*','LineWidth',2); % Dotted line
plot(list_PSNR,R_all_mca.R_all,'-.kx','LineWidth',2); % Dotted line
plot(list_PSNR,R_all_mmca.R_all,'--k^','LineWidth',2); % Dotted line
hold off
xlabel('PSNR in dB') 
ylabel('Correlation') 
legend('ICA','EFICA','JADE','MCA','MMCA','Location','southeast')
% --------------------------- %
figure;
hold on;
plot(list_PSNR,cri_all_ica.cri_all,'-kd','LineWidth',2);
plot(list_PSNR,cri_all_efica.cri_all,'--ko','LineWidth',2);
plot(list_PSNR,cri_all_jader.cri_all,':k*','LineWidth',2);
plot(list_PSNR,cri_all_mmca.cri_all,'-.k^','LineWidth',2); % Dotted line
hold off;
xlabel('PSNR in dB') 
ylabel('Mixing Matrix Criterion')
legend('ICA','EFICA','JADE','MMCA','Location','southeast')

