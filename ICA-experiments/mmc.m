% calculates the mixing matrix criterion
function [result] = mmc(A,piA,P)
% A: original mixing matrix
% piA: estimated demixing matrix
    [dim1,dim2] = size(A);
    dim = min(dim1,dim2);
    I = eye(dim,dim);
    pI = P*piA*A;
%     pI = abs(pI);
    for i = 1:dim
        pI(:,i) = pI(:,i)/sum(pI(:,i));
    end    
    result = norm(I - pI);
return
