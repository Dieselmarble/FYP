function SubBPtBP(A,k,l,Z)
%--------------------------------------------------------------------------
dim=size(A);  
m_subbp=floor(dim(1)/k);
n_subbp=floor(dim(2)/l);
%--------------------------------------------------------------------------
for i=1:m_subbp*n_subbp  %1:m_subbp*n_subbp
    row=ceil(i/n_subbp);
    column=n_subbp-(row*n_subbp-i);
    B=ones(k,l); 
    for j=1:l
        Zi=Z(i,:);
        B(:,j)=( Zi( (j-1)*k+1:k*j) )';
    end
    A(((row-1)*k+1):row*k,((column-1)*l+1):column*l)=B;
end
%--------------------------------------------------------------------------





    
    
    