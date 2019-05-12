%---???????????????????
function t=sidwt(t,wname)
if ~iscell(t{1,1})
    if ~iscell(t{1,2})
        if ~iscell(t{1,3})
            if ~iscell(t{1,4})
                [mi(1),ni(1)]=size(t{1,1});[mi(2),ni(2)]=size(t{1,2});
                [mi(3),ni(3)]=size(t{1,3});[mi(4),ni(4)]=size(t{1,4});
                m=min(mi);n=min(ni);
                t{1,1}=t{1,1}(1:m,1:n);
                t{1,2}=t{1,2}(1:m,1:n);
                t{1,3}=t{1,3}(1:m,1:n);
                t{1,4}=t{1,4}(1:m,1:n);
                t=idwt2(t{1,1},t{1,2},t{1,3},t{1,4},wname);
            else
                t{1,4}=sidwt(t{1,4},wname);
                [mi(1),ni(1)]=size(t{1,1});[mi(2),ni(2)]=size(t{1,2});
                [mi(3),ni(3)]=size(t{1,3});[mi(4),ni(4)]=size(t{1,4});
                m=min(mi);n=min(ni);
                t{1,1}=t{1,1}(1:m,1:n);
                t{1,2}=t{1,2}(1:m,1:n);
                t{1,3}=t{1,3}(1:m,1:n);
                t{1,4}=t{1,4}(1:m,1:n);
                t=idwt2(t{1,1},t{1,2},t{1,3},t{1,4},wname);
            end;
        else
            if iscell(t{1,4})
                t{1,4}=sidwt(t{1,4},wname);
            end;
            t{1,3}=sidwt(t{1,3},wname);
            [mi(1),ni(1)]=size(t{1,1});[mi(2),ni(2)]=size(t{1,2});
            [mi(3),ni(3)]=size(t{1,3});[mi(4),ni(4)]=size(t{1,4});
            m=min(mi);n=min(ni);
            t{1,1}=t{1,1}(1:m,1:n);
            t{1,2}=t{1,2}(1:m,1:n);
            t{1,3}=t{1,3}(1:m,1:n);
            t{1,4}=t{1,4}(1:m,1:n);
            t=idwt2(t{1,1},t{1,2},t{1,3},t{1,4},wname);
            
        end;
    else
        if iscell(t{1,4})
                t{1,4}=sidwt(t{1,4},wname);
        end;
        if iscell(t{1,3})
                t{1,3}=sidwt(t{1,3},wname);
        end;
        t{1,2}=sidwt(t{1,2},wname);
        [mi(1),ni(1)]=size(t{1,1});[mi(2),ni(2)]=size(t{1,2});
        [mi(3),ni(3)]=size(t{1,3});[mi(4),ni(4)]=size(t{1,4});
        m=min(mi);n=min(ni);
        t{1,1}=t{1,1}(1:m,1:n);
        t{1,2}=t{1,2}(1:m,1:n);
        t{1,3}=t{1,3}(1:m,1:n);
        t{1,4}=t{1,4}(1:m,1:n);
        t=idwt2(t{1,1},t{1,2},t{1,3},t{1,4},wname);
    end;
else
        if iscell(t{1,4})
                t{1,4}=sidwt(t{1,4},wname);
        end;
        if iscell(t{1,3})
                t{1,3}=sidwt(t{1,3},wname);
        end;
        if iscell(t{1,2})
                t{1,2}=sidwt(t{1,2},wname);
        end;
        t{1,1}=sidwt(t{1,1},wname);
        [mi(1),ni(1)]=size(t{1,1});[mi(2),ni(2)]=size(t{1,2});
        [mi(3),ni(3)]=size(t{1,3});[mi(4),ni(4)]=size(t{1,4});
        m=min(mi);n=min(ni);
        t{1,1}=t{1,1}(1:m,1:n);
        t{1,2}=t{1,2}(1:m,1:n);
        t{1,3}=t{1,3}(1:m,1:n);
        t{1,4}=t{1,4}(1:m,1:n);
        t=idwt2(t{1,1},t{1,2},t{1,3},t{1,4},wname);
end;
return;