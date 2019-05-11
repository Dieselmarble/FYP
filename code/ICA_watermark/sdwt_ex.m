%---??????????????
function  t=sdwt_ex(node,wname,tkey)
[A,H,V,D]=dwt2(node,wname);
AHVD={A,H,V,D};
for i=1:length(tkey)
    if iscell(tkey{1,i})
        t{1,i}=sdwt_ex(AHVD{1,i},wname,tkey{1,i});
    else
        if tkey{1,i}~=-1
            t{1,i}=node;
        else
            t{1,i}={};
        end;
    end;
end;