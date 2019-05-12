%Short demo
%


%%%%%%%%%%%
%Simple use
%%%%%%%%%%%
n=5; %number of sources
T=10000; %Number of samples
s=rand(n,T);  %sources
A=randn(n,n); %mixing matrix
x=A*s;  %mixtures
B=ogwe(x);   %Separation-projection matrix. 
             %With no more parameter than x it
             %executes the SICA method with no 
             %dimension reduction

disp('Product of Separation and Mixing matrix')
B*A  %should be a non-mixing matrix

%%%%%%%%%%%%%%%%%%%
%Other posibilities
%%%%%%%%%%%%%%%%%%%
%A) reduction to m dimension by space projection
m=4;
B=ogwe(x,m);  

%B) Selecting algorithm
wxi=3/7; wkur=0; %SICA method
wxi=0; wkur=+1;  %ML, EML, MK, SKSE methods for positive kurtosis
                 %Also wkur=-1; for negative kurtosis
wxi=1/3; wkur=0; %AML method
wxi=0; wkur=0;   %AEML method
wxi=1/2; wkur=0; %MASSFOC method
               
B=ogwe(x,wkur,wxi);  

%C) Selecting algorithm with dimension reduction
B=ogwe(x,m,wkur,wxi);  