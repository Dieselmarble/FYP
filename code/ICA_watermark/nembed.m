%---???????????
function [wnode,nodekey]=nembed(nd,b)
[m1,n1]=size(nd);
[m2,n2]=size(b);
node=nd;
node1=node(1:m2*n2);
node2=node(m2*n2+1:m1*n1);
b=b(1:m2*n2);

%---????---------------------------------------------------------------
S=[];
for i=1:2
    switch i
    case 1, s=double(node1);
    case 2, s=double(b);
    end
    s=s-mean(s);   % ??
    s=s/std(s,1);  % ???
    S=[S; s];
end
A=rand_orth(2);
% A=rand(2,2);
X=A*S;
x1=X(1,:);
x1=[x1,node2];
x2=X(2,:)
nodekey=reshape(x2,m2,n2);
%--------------------------------------------------------------------------
%---??-------------------------------------------------------------------
% x1=node1+b*0.02;
% x1=[x1 node2];
% nodekey=[0.02,-1];
%--------------------------------------------------------------------------
wnode=reshape(x1,m1,n1);

return;