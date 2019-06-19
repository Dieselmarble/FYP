% calculates the mixing matrix criterion
% A: original mixing matrix
% piA: estimated demixing matrix
function [result] = mmc(A,piA,P)
    [dim1,dim2] = size(A);
    dim = min(dim1,dim2);
    I = eye(dim,dim);
    pI = P*piA*A;
    for i = 1:dim
        pI(:,i) = abs(pI(:,i))/max(abs(pI(:,i))); % scale the difference
    end    
%     pI = abs(pI);
%     for i = 1:dim
%         pI(:,i) = pI(:,i)/sum(pI(:,i));
%     end    
%     result = sum(sum(pI-I)/(dim*(dim-1));
    result = norm(abs(I-pI),'fro');
return
