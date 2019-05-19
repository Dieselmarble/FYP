function index=max_cov(Q)
    dim=size(Q);
    a=cov(Q(1,:));
    index=1;
    for i=1:dim(1)
        if cov(Q(i,:))<a 
            a=cov(Q(i,:));
            index=i;
        end
end
