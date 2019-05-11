%---????????????????----
function  [t,tkey]=sdwt(node,wname,m2,n2,e0)
[A,H,V,D]=dwt2(node,wname);
[mi,ni]=size(A);
if mi*ni<m2*n2
    t=node;tkey=1;
    return;
end;
ea=(sum(sum(A.^2)))/(mi*ni);
eh=(sum(sum(H.^2)))/(mi*ni);
ev=(sum(sum(V.^2)))/(mi*ni);
ed=(sum(sum(D.^2)))/(mi*ni);
if ea<e0&ea>0.1*e0
    [A,ka]=sdwt(A,wname,m2,n2,e0);
else
    ka=-1;
end;
if eh<e0&eh>0.1*e0
    [H,kh]=sdwt(H,wname,m2,n2,e0);
else
    kh=-1;
end;
if ev>e0&ev>0.1*e0
    [V,kv]=sdwt(V,wname,m2,n2,e0);
else
    kv=-1;
end;
if ed<e0&ed>0.1*e0
    [D,kd]=sdwt(D,wname,m2,n2,e0);
else
    kd=-1;
end;
    t={A,H,V,D};
    tkey={ka,kh,kv,kd};

return;
