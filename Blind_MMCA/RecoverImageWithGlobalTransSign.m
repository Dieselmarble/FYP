function [yout yout1]=RecoverImageWithGlobalTransSign(y,x,mu,lambda2,D,CoefMatrix,blocknum,stepsize,dc,transform)
% ========================================================
% ========================================================
N=size(y,1); 
n=sqrt(size(D,1)); 
yout=zeros(size(y)); 
Weight=zeros(size(y));

% a kind of soft thresholding
switch transform
    case 'dct'
        newterm = idct2(sign(dct2(x)));
    case 'wavelet'
        qmf = MakeONFilter('Battle',5);
        newterm = IWT2_PO(sign(FWT2_PO(x,1,qmf)),1,qmf);
end
i=1; j=1;
data= D*CoefMatrix;
if dc ~= 0
    data = add_dc(data,dc,'columns');
end
for k=1:1:blocknum,
    patch=reshape(data(:,k),[n,n]); 
    yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch; 
    Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1; 
    if i<N-n+1 
        i=i+stepsize(1); 
    else
        i=1; j=j+stepsize(2); 
    end;
end;
yout1 = yout;
yout=(mu*yout+y-lambda2*newterm)./(mu*Weight+1); 
return;
