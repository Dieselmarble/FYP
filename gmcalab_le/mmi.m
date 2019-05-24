% calculates the mixing matrix criterion
function [result] = mmi(A,piA,P)
    [dim1,dim2] = size(A);
    dim = min(dim1,dim2);
    result = eye(dim,dim);
    result = result - P*piA*A;
end
