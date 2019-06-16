function [yout yout1 Weight]=Recover1DSignalWithGlobalTrans2(y,mu,D,CoefMatrix,blocknum,stepsize,dc)
% ========================================================
% ========================================================
N=size(y,2); 
n=(size(D,1)); 
yout=zeros(size(y)); 
Weight=zeros(size(y));
i=1; j=1;
data= D*CoefMatrix;
if dc ~= 0
    data = add_dc(data,dc,'columns');
end
for k=1:1:blocknum,
    seg=data(:,k)'; 
    yout(i:i+n-1)=yout(i:i+n-1)+seg; 
    Weight(i:i+n-1)=Weight(i:i+n-1)+1; 
    if i<N-n+1 
        i=i+stepsize; 
    else
        i=1;
    end;
end;
yout1 = yout;
yout=(mu*yout+y); 
return;
