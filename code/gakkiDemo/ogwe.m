function B = ogwe(x,m,wkur,wxi)
%
%B =  ogwe(x,m,wkur,wxi)
%
%B=OGWE(x) is the optimized SICA BSS/ICA real algorithm 
%   returning B for mixtures in the rows of x using the
%   SICA with parameters. Thus, estimated sources 
%   (independent components) are y=B*x
%
%B=OGWE(x,m)
%   Parameter m allows a reduction on the dimension 
%   via subspace projection.
%   
%B=OGWE(x,wkur,wki) is the Optimized Genaralized Weigthed 
%   estimator. Some particular cases of the GWE estimator 
%   are:
%    wxi=3/7, wkur=0   -->   SICA method
%    wxi=1/2, wkur=0   -->   MASSFOC method
%    wxi=1/3, wkur=0   -->   AML method
%    wxi=0, wkur=0     -->   AEML method
%    wxi=0, wkur=+-1   -->   ML, EML, MK, SKSE methods
%                            wkur= 1 for positive kurtosis sources and 
%                            wkur=-1 if negative kurtosis sources
%B=OGWE(x,m,wkur,wki)
%   Parameter m allows a reduction on the dimension 
%   via subspace projection.
%
%  You may refer to
%  a) Vicente Zarzoso, Juan J. Murillo-Fuentes, Rafael Boloix-Tortosa, Asoke K. Nandi. 
%     “Optimal pairwise fourth-order independent component analisis”. 
%     IEEE Trans. On Signal Processing, vol. 54, pp. 3049-3063. 2006.
%  b) Juan J. Murillo-Fuentes, Francisco J. González-Serrano. 
%     “A sinusoidal contrast function for the blind separation of 
%     statistically independent sources”. IEEE Trans. On Signal Processing, 
%     vol. 52, pp. 3459-3463. 2004.
%  for further details on the algorithm.
%
%  Please report any problem or bug to murillo@us.es
%
% Version 1.0 21/Jan/2003
% Copyright (c): Juan Jose Murillo-Fuentes and Francisco Javier Gonzalez-Serrano.
%
%Permission is granted for anyone to copy, use, or modify these programs for
%purposes of research or education, provided this copyright notice is retained,
%and note is made of any changes that have been made.
%
%These programs are distributed without any warranty, express or
%implied. As these programs were written for research purposes only, they
%have not been tested to the degree that would be advisable in any
%important application.  All use of these programs is entirely at the
%user's own risk.
%
%

   

[n,T]=size(x);
verbose=0 ;	% Set to 0 for quiet operation

switch(nargin)
case 1 
    wxi=3/7;
    wkur=0;
    m=n;
case 2
    if m-floor(m)>0 | m<=0 fprintf('m should be a positive integer!\n'), return,end    
    if m>n ,    fprintf('No more sources than sensors here!\n'), return,end    
    if verbose, fprintf('Looking for %d sources\n',m); end ;
    wxi=3/7;
    wkur=0;
case 3    
    wxi=m;
    wkur=wxi;
    m=n;
case 4    
    if m-floor(m)>0 | m<=0 fprintf('m should be a positive integer!\n'), return,end        
    if m>n ,    fprintf('No more sources than sensors here!\n'), return,end    
    if verbose, fprintf('Looking for %d sources\n',m); end ;
end



if verbose, fprintf('Zero mean\n'); end 
x	= x - mean(x')' * ones(1,T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% whitening 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose, fprintf('Whitening\n'); end
[V,D]= eig((x*x')/T); 
[vals,idxs]	= sort(diag(D));
mmvar = n-m+1:n; 
V=V(:,idxs(mmvar))';
D=D(idxs(mmvar),idxs(mmvar));
W = sqrt(diag(1./diag(D)))  * V;
x = W*x;  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Contrast minimization to compute unitary Matrix R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose, fprintf('Contrast optimization: Unitary matrix computation\n'); end

if m==2,K=1;else,K=1+round(sqrt(m));end;  % K= max number of sweep  Comon
seuil=pi/360; %A fixed allowed angle error as threshold or,...
%seuil	= 1/sqrt(T)/100; % A statistically significant threshold

burden=1;
sweep = 0;
GivensRotations=0;
R	= eye(m) ; % the initial rotation matrix


%%%%%%%%%%CHOSE between Not Initialized SICA or Initialized SICA

a=123.75;b=-237.5;
y=a*n+b;

%fprintf('sources and samples, n=%d and T=%d, flops %d\n',n,T,y)

if n<15 & T>y
    INITIALIZATION=1; % 'ISICA';  % disp('ISICA')
else 
    INITIALIZATION=0; %'SICA';   % disp('SICA')
end



if INITIALIZATION  %%%%%%%%%%%%%%ISICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nm=m*(m+1)/2;  %This is the order of the Momment Matrix 
    X4=zeros(nm,nm);  %Initialitation of the Momment Matrix 
    ix1=1;   %Indexes to fill the entries of the Momment Matrix
    ix2=m; 
    
    %The Momment Matrix is used to compute some momments of order 4. When a new givens angle is 
    %computed, we do not rotate the data by this angle and recompute the momment but we multiply
    %the momment matrix by the rotation matrix and compute the momment. As an example, if we are
    %to calculate Mpppp, then we write Mppp=Rpp*X4*Rpp'. Rpp is a vector containing the non-repiting
    %entries of the product R(p,:)'*R(p,:). We initialize some vectors first
    
    %%%%%%%%%%%%%%%%%%
    %Rotation vectors and Moment Matrix initialization
    
    if verbose, fprintf(' Initializing moment matrix and rotation vectors\n'); end
    Rpp=zeros(1,nm); 
    Rqq=zeros(1,nm);
    Rpq=zeros(1,nm);
    %Xpq1=zeros(1,T);
    Xpq=zeros(1,T);
    
    p1=[];
    q1=[];
    ix=1;
    for p_1=1:m
        sccero=.5;
        scuno=1;
        for q_1=p_1:m
            sc1(ix)=scuno;
            sc0(ix)=sccero;
            p1=[p1 p_1];
            q1=[q1 q_1];
            scuno=2;
            sccero=1;
            ix=ix+1;
        end
    end
    
    %Number of blocks
    Nb=m;
    dim=sum(m:-1:1);
    ixq=[];vq=1;  
    ixk1=[];
    ixr=[];
    for l=Nb:-1:1
        ixq=[ixq vq*ones(1,l)];
        ixk1=[ixk1 [vq:Nb]];
        vq=vq+1;   
        ixr=[ixr sum(Nb:-1:l+1)+1];
    end
    ixr=[ixr ixr(m)];
    
    Xqk1=zeros(ixr(m),T);
    Xqk1=x(ixq,:).*x(ixk1,:);
    for p=1:Nb
        k2=p:m;
        idc=ixr(p)+k2-p;
        idr=ixr(k2(1)):dim;
        X4(idr,idc)=Xqk1(idr,:)*(Xqk1(idc,:))'/T;      
        X4(idc,idr)=X4(idr,idc)';
        idc=ixr(p);
        for k2=p:m    
            for k3=idc:ixr(k2)-1
                ixcoj=[sort([ixk1(k3) k2])];
                col=ixr(p)-p+ixq(k3);
                row=ixr(ixcoj(1))-ixcoj(1)+ixcoj(2);
                X4(k3,idc)=X4(row,col);
                X4(idc,k3)=X4(k3,idc);
            end
            idc=idc+1;   
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%
    %Contrast minimization
    
    %while burden, burden=0;
    %while sweep < K
    while burden & sweep < K, burden=0;
        
        if verbose, fprintf(' Sweep #%d\n',sweep); end
        sweep=sweep+1;
        
        for p=1:m-1,   
            for q=p+1:m,
                
                %We compute here the rotations vectors needed to calculate the momments
                Rpp= sc1.*R(p,p1).*R(p,q1); %               
                Rqq= sc1.*R(q,p1).*R(q,q1); %               
                Rpq= sc0.*(R(p,p1).*R(q,q1)+R(p,q1).*R(q,p1));
                
                %Then we compute the momments
                RppX4=Rpp*X4;
                X4Rqq=X4*Rqq';
                Mppqq=RppX4*Rqq';
                Mpppp=RppX4*Rpp';
                Mqqqq=Rqq*X4Rqq;
                Mpqqq=Rpq*X4Rqq;
                Mpppq=RppX4*Rpq';         
                
                r4=Mpppp+Mqqqq+2*Mppqq;  %r4=mean(r4);         %r4=Mpppp+(Rqq+2*Rpp)*X4Rqq;
                sxpxqr2=2*(Mpqqq+Mpppq);      %  sx1x2r2=2*mean(x1x2.*r2); %%%%
                
                r4c4phi=r4-8*Mppqq;
                r4s4phi=8*Mpppq-2*sxpxqr2;
                r4c2phi=Mpppp-Mqqqq;  %r4c2phi=2*(Mpppp+Mppqq)-r4;
                r4s2phi=sxpxqr2;
                
                %%% computation of Givens angle
                
                r4c2phi2=r4c2phi^2;
                r4s2phi2=r4s2phi^2;
                
                r4_c=r4-8;  %c=1.5/.1875
                
                xi4=r4c4phi+j*r4s4phi;
                xi2_2=(r4c2phi+j*r4s2phi)^2;
                theta=angle(wkur*xi4 + (1-abs(wkur))*(wxi*r4_c*xi4+(1-wxi)*xi2_2) )/4; %GWE                
                
                %%% Givens update
                GivensRotations=GivensRotations+1;
                if abs(theta) > seuil, burden=1 ;
%                    updates = updates + 1;
                    c	= cos(theta); 
                    s	= sin(theta);
                    G	= [ c s ; -s c ] ;
                    pair = [p;q] ;
                    R(pair,:)=G*R(pair,:);
                end
            end%%of the loop on q
        end%%of the loop on p
    end
    
else  %%%%%%%%%%%%%%NISICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %while burden, burden=0;
    %while sweep < K
    while burden & sweep < K, burden=0;
        
        if verbose, fprintf(' Sweep #%d\n',sweep); end
        sweep=sweep+1;
        
        for p=1:m-1,
            for q=p+1:m,
                Mpp=x(p,:).*x(p,:);
                Mqq=x(q,:).*x(q,:);
                Mpq=x(p,:).*x(q,:);
                Mpppp=Mpp*Mpp'/T;          
                Mqqqq=Mqq*Mqq'/T;
                Mpqqq=Mpq*Mqq'/T;
                Mpppq=Mpp*Mpq'/T;
                Mppqq=Mpp*Mqq'/T;
                
                r4=Mpppp+Mqqqq+2*Mppqq;  
                sxpxqr2=2*(Mpqqq+Mpppq); 
                
                r4c4phi=r4-8*Mppqq;
                r4s4phi=8*Mpppq-2*sxpxqr2;
                r4c2phi=Mpppp-Mqqqq;  %r4c2phi=2*(Mpppp+Mppqq)-r4;
                r4s2phi=sxpxqr2;
                
                %%% computation of Givens angle
                
                r4c2phi2=r4c2phi^2;
                r4s2phi2=r4s2phi^2;
                
                r4_c=r4-8;  
                xi4=r4c4phi+j*r4s4phi;
                xi2_2=(r4c2phi+j*r4s2phi)^2;
                theta=angle(wkur*xi4 + (1-abs(wkur))*(wxi*r4_c*xi4+(1-wxi)*xi2_2) )/4; %GWE                
                
                
                %%% Givens update
                GivensRotations=GivensRotations+1;                
                if abs(theta) > seuil, burden=1 ; 
                    c	= cos(theta); 
                    s	= sin(theta);
                    G	= [ c s ; -s c ] ;
                    pair = [p;q] ;
                    R(pair,:)=G*R(pair,:);
                    x([p q],:)=G*x([p q],:);
                end%%of the if
            end%%of the loop on q
        end%%of the loop on p
    end
end
    if verbose, fprintf('%d givens angle computed\n',GivensRotations); end    


    
    
%%%%%%%%%%%%%%%%%%%
%Final computations
%%%%%%%%%%%%%%%%%%%

    %Final result is whitening plus unitary transformation
    B	= R* W ;
   
    % We rearrange matrix B to get components with higher variance first
    if verbose, fprintf('Sorting components\n'); end
    P=inv(D)*R';
    [vals,idxs]	= sort(sum(P.*P)) ;
    B		= B(idxs,:);

    
return ;
