%this function is refer to https://ieeexplore.ieee.org/document/941340
function Z = DvdBptSubBp(A,k,l) %k,l is dimension of the watermark
dim=size(A); %dimension 
m_subbp=floor(dim(1)/k); %p on paper
n_subbp=floor(dim(2)/l); %q on paper
Z=ones(m_subbp*n_subbp,k*l); %C_pq on paper; dimension of Z is equal to dimension of A
%------------------------------------------
for i=1:m_subbp*n_subbp
    row=ceil(i/n_subbp); % row = m_subbp; This is equivalent to two loops; ceil is opposite of 'floor()'
    column=n_subbp-(row*n_subbp-i);
    temp=A((row-1)*k+1: row*k, (column-1)*l+1: column*l); %k+1:2k ;l+1:2l
    Z(i,:)= temp(:)';
end
