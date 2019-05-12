%--???????????
function [wi,map]=nextract(tnd,tkey)
[m2,n2]=size(tkey);
x1=tnd(1:m2*n2);
x2=tkey(1:m2*n2);
X=[x1;x2];
B=jadeR(double(X),2);
Y=B*double(X);
y2=Y(2,:);
wi=y2;
map=[m2,n2];
return;