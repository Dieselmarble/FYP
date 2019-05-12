%---?????????????
function [t,tkey]=embed(t,tkey,b)
[m2,n2]=size(b);
for i=1:4
    if ~iscell(tkey{1,i})
        if tkey{1,i}==1
        [t{1,i},tkey{1,i}]=nembed(t{1,i},b);
        end;        
    else
        [t{1,i},tkey{1,i}]=embed(t{1,i},tkey{1,i},b);
    end; 
end;
return;