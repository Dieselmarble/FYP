function [yout yout1 Weight]=RecoverImageWithGlobalTrans2(y,mu,D,CoefMatrix,blocknum,stepsize,dc)
% ========================================================
% ========================================================
N=size(y,1); 
n=sqrt(size(D,1)); 
yout=zeros(size(y)); 
Weight=zeros(size(y));
i=1; j=1;a=0;k=1;
data= D*CoefMatrix;
if dc ~= 0
    data = add_dc(data,dc,'columns');
end
while k<=blocknum,
    patch=reshape(data(:,k),[n,n]); 
    if i+n-1 <= size(yout,1) && j+n-1 <= size(yout,2)
        yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch; 
        Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1; 
    else
        a=a+1;k=k-1;
    end
    if i<N-n+1 %&& i+stepsize(2)+n-1 < size(yout,2)
        i=i+stepsize(2);  
    else
        i=1; 
        %if j+stepsize(1)+n-1 < size(yout,1)
            j=j+stepsize(1); 
%        end
    end;
    k=k+1;
end;
yout1 = yout;
yout=(mu*yout+y); 
return;
