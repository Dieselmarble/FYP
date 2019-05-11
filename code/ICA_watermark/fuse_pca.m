%---PCA?????????
function Y = fuse_pca(M1, M2)
[z1 s1] = size(M1);
[z2 s2] = size(M2);
if (z1 ~= z2) | (s1 ~= s2)
  error('Input images are not of same size');
end;
[V, D] = eig(cov([M1(:) M2(:)]));
if (D(1,1) > D(2,2))
  a = V(:,1)./sum(V(:,1));
else  
  a = V(:,2)./sum(V(:,2));
end;
Y = a(1)*M1+a(2)*M2;