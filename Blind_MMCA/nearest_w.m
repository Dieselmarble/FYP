function [A,M]=nearest_w(Wnew,Wold);
n=size(Wnew,2);
Sign=zeros(n,n);
for i=1:n
     for j=1:n
         [W(i,j) Sign(i,j)]=min([sum((Wnew(:,i)-Wold(:,j)).^2) sum((Wnew(:,i)+Wold(:,j)).^2)]);
    end
end
Sign=-sign(Sign-1.5);
W=max(max(W))-W;
M=maxmatching(W);
M=M.*Sign;
A=Wnew*M; 