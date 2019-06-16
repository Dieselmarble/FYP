function [X C] = simult_sparse_coding(D, Y, d, k, find_X)
L = size(Y,2); [N,K] = size(D);
n = sqrt(sum(D.*D));
Dn =ones(N,1)*(1./n).*D;
dd = false(max(d), K);
for i = 1:max(d), 
    dd(i,d==i) = true; % group by blocks
end;
[X C] = simult_BMMP(Dn, Y, dd, k, find_X); % batch matching pursuit
if find_X, 
    X = (1./n')*ones(1,L).*X; 
end
