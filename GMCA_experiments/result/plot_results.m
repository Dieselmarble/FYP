close all ; clear all; clc;
% ---------------------------- %
cri_all_gmca = load('cri_gmca.mat');
cri_all_efica = load('cri_all_efica.mat');
cri_all_jader = load('cri_all_jader.mat');
cri_all_mmca = load('cri_all_MMCA.mat');

% ---------------------------- %
R_all_gmca = load('R_gmca.mat');
R_all_efica = load('R_all_efica.mat');
R_all_jader = load('R_all_jader.mat');
R_all_mmca = load('R_all_MMCA.mat');
% ---------------------------- %
list_PSNR = [3:5:40];
figure;
hold on
plot(list_PSNR,R_all_gmca.R_all,'-kd','LineWidth',2); % Solid line (default)
plot(list_PSNR,R_all_efica.R_all,'--ko','LineWidth',2); % Dashed line
plot(list_PSNR,R_all_jader.R_all,':k*','LineWidth',2); % Dotted line
plot(list_PSNR,R_all_mmca.R_all,'-.kx','LineWidth',2); % Dotted line
hold off
xlabel('PSNR in dB') 
ylabel('Correlation') 
legend('GMCA','EFICA','JADE','MMCA','Location','southeast')
% --------------------------- %
figure;
hold on;
plot(list_PSNR,cri_all_gmca.cri_all,'-kd','LineWidth',2);
plot(list_PSNR,cri_all_efica.cri_all,'--ko','LineWidth',2);
plot(list_PSNR,cri_all_jader.cri_all,':k*','LineWidth',2);
plot(list_PSNR,cri_all_mmca.cri_all,'-.k^','LineWidth',2); % Dotted line
hold off;
xlabel('PSNR in dB') 
ylabel('Mixing Matrix Criterion')
legend('GMCA','EFICA','JADE','MMCA','Location','southeast')

