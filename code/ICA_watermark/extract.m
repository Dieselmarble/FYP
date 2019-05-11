%---????????
function [w,map]=extract(t,tkey);
w=[];
for i=1:4
    if ~iscell(tkey{1,i})
        if tkey{1,i}~=-1
            [wi,map]=nextract(t{1,i},tkey{1,i});
        else
            wi=[];
        end; 
    else
        [wi,map]=extract(t{1,i},tkey{1,i});
    end; 
    if length(wi)>0
        w=[w;wi];
    end;
end;
return;