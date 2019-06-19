function A = Optisens(D, M)
%Sensing matrix design for sparse coding

[N K] = size(D);
[U,S,V] = pca(D*D',M);
A = sqrt(inv(S))*U';


