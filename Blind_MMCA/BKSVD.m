function [D2,X2,d2] = BKSVD(Y,K,k,s,d0,max_it,D1);
%Run BKSVD+SAC

% initialize dictionary with random columns of norm 1
% D1 = randn(size(Y,1),K);
D1 = D1 ./  repmat(sqrt(sum(D1.^2)),size(Y,1),1);
D2 = D1;
X1 = []; C1 = [];
% h=waitbar(0,'BK-SVD & SAC');
for i = 1:max_it
    %ks-sparse reprentations of Y w.r.t D2 are calculated
    %(0 means only C2s is calculated, but not X2s, this is all we need for SAC)
    %C2s of size LxK contains true at (i,j) if X0(j,i)!=0
    [X2s C2s] = simult_sparse_coding(D2,Y,d0,k*s,0);
    %SAC-finds block structure d2 given C2s and max block size s
    d2 = sparse_agg_clustering(C2s, s);
    %k-block-sparse reprentations X2 over d2 of Y w.r.t D2 are calculated
    %C2 of size LxB contains true at (i,j) if X2(:,i)  'uses' block nr j
    [X2 C2] = simult_sparse_coding(D2,Y,d2,k,1);
    %KSVD - updates every atom in D2 and nonzero values in X2 to minimize representation error
    [X2 D2] = KSVD(X2, D2, Y, d2, C2);
%     waitbar(i/max_it);
end
% close(h);
end