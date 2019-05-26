% calculates the mixing matrix criterion
function [result] = mmc(A,piA,P)
    [dim1,dim2] = size(A);
    dim = min(dim1,dim2);
    I = eye(dim,dim);
    pI = P*piA*A;
    for i = 1:4
        pI(:,i) = pI(:,i)/sum(pI(:,i));
    end    
    result =  norm(I - pI);
end
