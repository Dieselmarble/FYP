function Ahat = colnorm(B,n)

% A = colnorm(B,n);
% If only one argument is used the columns of matrix B is normalized
% Otherwise, a column normed matrix B*n is created

if nargin > 1
    B = randn(B,n);
end     

% normalizing the columns of Dictionary
for i = 1:size(B,2)
    B(:,i) = B(:,i)/norm(B(:,i));
end
Ahat = B;
end        
