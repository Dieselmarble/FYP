function Z = DvdBptSubBp(A,k,l)
dim=size(A); %dimension 
m_subbp=floor(dim(1)/k);
n_subbp=floor(dim(2)/l);
%------------------------------------------
Z=ones(m_subbp*n_subbp,k*l);
%------------------------------------------
for i=1:m_subbp*n_subbp
    row=ceil(i/n_subbp);
    column=n_subbp-(row*n_subbp-i);
    temp=A(((row-1)*k+1):row*k,((column-1)*l+1):column*l);
    Z(i,:)=(temp(:))';
end
