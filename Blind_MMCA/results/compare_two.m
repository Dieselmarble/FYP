close all
ksvd = [5.9319e+03,4.0494e+03,3.4322e+03,3.1677e+03,3.0284e+03,2.7654e+03;
     3.9792e+03,3.3865e3,3.0647e3,2.7270e3, 2.4699e+03,2.1415e+03;
    3.3636e+03,2.9656e3,2.5638e3,2.5475e3, 2.0406e+03,2.1898e+03;
    3.0900e+03,2.6525e3,2.6053e3,2.4809e3, 2.0940e+03,1.9050e+03;
    2.8856e+03,2.4214e+03,2.2733e+03,2.0863e+03,1.9877e+03,1.8098e+03
    2.7369e+03,2.3378e+03,2.1791e+03,2.0571e+03,1.8150e+03,1.5425e+03];

bksvd=[5.9319e+03,5.5859e+03,5.7371e+03,5.4835e+03,5.3238e+03,5.1695e+03;
     3.9684e+03,3.4703e3,3.1834e3,3.0298e3, 2.8209e+03,2.6118e+03;
    3.4143e+03,2.8892e3,2.6569e3,2.4676e3, 2.3138e+03,2.2513e+03;
    3.0910e+03,2.5149e3,2.3191e3,2.1887e3,2.0197e+03,2.0841e+03;
    2.9431e+03,2.1501e+03,2.0081e+03,1.8471e+03,1.9193e+03,1.8324e+03;
    2.7370e+03,1.9430e+03,1.7546e+03,1.9098e+03,1.6645e+03,1.6680e+03];

ax1 = 1:1:6;
ax2 = 1:1:6;
% s = contour(ax1,ax2,ksvd,'ShowText','on')
s = heatmap(ksvd,'ColorLimits',[0 3500]);
% set(gca,'XTick',1:1:4)
% set(gca,'YTick',1:1:4)
s.xlabel('Maximal Block Size')
s.ylabel('Block Sparsity Level')
colormap gray
% colorbar
figure
% s = contour(ax1,ax2,bksvd,'ShowText','on')
s = heatmap(bksvd,'ColorLimits',[0 3500]);
% set(gca,'XTick',1:1:4)
% set(gca,'YTick',1:1:4)
s.xlabel('Maximal Block Size')
s.ylabel('Block Sparsity Level')
colormap gray
% colorbar

% 1.3334e+03
% 
% 
% 1.4862e+03
