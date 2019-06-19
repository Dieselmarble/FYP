close all ; clear all; clc;
% ---------------------------- %
error_all_ksvd = load('errors_ksvd.mat');
error_all_bksvd = load('errors_bksvd.mat');
error_all_gmca = load('error_gmca.mat');

% ---------------------------- %
list_iter = [1:100];
figure;
hold on
plot(list_iter,error_all_ksvd.his1,'-k','LineWidth',2); % Solid line (default)
plot(list_iter,error_all_bksvd.his1,'--k','LineWidth',2); % Dashed line
plot(list_iter,error_all_gmca.his1,':k','LineWidth',2); % Dotted line
hold off
xlabel('Iteration Number') 
ylabel('Total Mean Square Error') 
legend('K-SVD','BK-SVD','Location','southeast')
% --------------------------- %

