%% Init 
clear all;
% close all;
%clc;

%% First we need to train a dictionary

%% patch parameters
PATCH.ROWS = 8;      % number of rows in each patch
PATCH.COLS = 8;      % number of columns in each patch
PATCH.COLORS = 1;    % use gray-scale or RGB
PATCH.SZ = PATCH.ROWS*PATCH.COLS*PATCH.COLORS;
N = PATCH.ROWS*PATCH.COLS;      % number of pixels in patch


%% load an image for training
ima = imread('barbara256.png');
if size(ima,3)>1, ima = rgb2gray(ima); end


%% extract all distinct patches in training image
Y = im2col(double(ima), [PATCH.ROWS PATCH.COLS], 'distinct');

%% Create Dictionary with training signals: 

N = PATCH.ROWS*PATCH.COLS;      % number of pixels in patch
s = 6; %dimension of blocks
k = 7;
K = 96;           %Number of columns in dictionary
L = size(Y,2);    %Number of training signals. should be thousands. Too many will be too heavy for PC
max_it = 100;       %Nr of iterations for the algorithm to converge. 
B = K/s;          %Number of blocks in D. Should be an order of magnitude higher than the 
                  %block sparsity k, otherwise the representations wont be sparse
d0 = 1:K;         %block structure with K blocks of size 1 (i.e. block structure ignored)
d = repmat(1:B, s,1); d = d(:)';  %block structure with B blocks of size s
M = 16; %number of samples after compression

%% train dictionaries
% D0 = randn(size(Y,1),K);
% D0 = D0 ./  repmat(sqrt(sum(D0.^2)),size(Y,1),1);
% Dictionary = D0;
% Data=Y;
% patchsz=8;
% sigma=0;
% const=sqrt(1.15);
% for i = 1:max_it
%     % Batch-OMP
%     CoefMatrix=omp2(Dictionary'*Data,sum(Data.*Data),Dictionary'*Dictionary,patchsz*sigma*const);
%     for j=1:1:size(Dictionary,2) % implement the K-SVD algorithm
%         relevantDataIndices=find(CoefMatrix(j,:)); %Find indices of nonzero elements.
%         if ~isempty(relevantDataIndices)
%             tmpCoefMatrix=CoefMatrix(:,relevantDataIndices);
%             tmpCoefMatrix(j,:)=0;
%             errors=Data(:,relevantDataIndices)-Dictionary*tmpCoefMatrix;
%             [betterDictionaryElement,singularValue,betaVector]=svds(errors,1);
%             CoefMatrix(j,relevantDataIndices)=singularValue*betaVector';
%             Dictionary(:,j)=betterDictionaryElement;
%         end;
%     end;
% end
% Ycomp2 = Dictionary*CoefMatrix;
% YcompImage2 = col2im(Ycomp2, [PATCH.ROWS PATCH.COLS], [size(ima,1) size(ima,2)], 'distinct');
% imshow(YcompImage2,[]);
% return
% % -------------------------------------------------- %
% initialize dictionary with random columns of norm 1
D0 = randn(size(Y,1),K);
D0 = D0 ./  repmat(sqrt(sum(D0.^2)),size(Y,1),1);
D1 = D0;
X1 = []; C1 = [];
for i = 1:max_it
    %block sparse coding
    %k*s-sparse reprentations of Y w.r.t D1 are calculated
    %(d0= [1,2,...,K] means sparse and not block sparse)
    %1 means X1 is explicitly calculated
    [X1 C1] = simult_sparse_coding(D1,Y,d0,k*s,1);
    %KSVD - updates every atom in D1 and nonzero values in X1 to minimize representation error
    [X1 D1] = KSVD(X1, D1, Y, d0, C1);
    % KSVD method and BKSVD method use same number of iterations:
%     if (i==max_it)
%         D2=D1; % save KSVD results after half the iterations
%     end
end
e1 = norm(Y-D1*X1,'fro')

%  D1 - result of KSVD
%  D2 - result of BKSVD+SAC
%  d2 - recovered block structure for BKSVD+SAC
%  D3 - result of BKSVD without SAC
[D2,X2] = BKSVD_old(Y,K,k,s,d0,max_it,D0);
e2 = norm(Y-D2*X2,'fro') %BKSVD error
% Dictionary = rand(64,96);
% [D2,CoefMatrix,e2] = BKSVD(Y,K,k,s,d0,max_it,Dictionary);
% dictimg = showdict(D2,[1 1]*8,round(sqrt(256)),round(sqrt(256)),'lines','highcontrast');
% imshow(imresize(dictimg,2,'nearest'));
% title('Trained dictionary');
Ycomp2 = D2*X2;
YcompImage2 = col2im(Ycomp2, [PATCH.ROWS PATCH.COLS], [size(ima,1) size(ima,2)], 'distinct');
imshow(YcompImage2,[]);
return;
%% learn A using Sapiro's algorithm
A = Optisens(randn(N,K),M);


%% load an image for testing
ima = imread('cute_baby.jpg');
if size(ima,3)>1, ima = rgb2gray(ima); end

%%%%%%%%%%%%%%%%%%%%% Distinct patches %%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract all distinct patches in training image
Y = im2col(double(ima), [PATCH.ROWS PATCH.COLS], 'distinct');

%% reconstruct signals from compressed vectors
X1 = simult_sparse_coding(A*D1,A*Y,d0,k*s,1);
X2 = simult_sparse_coding(A*D2,A*Y,d2,k,1);
X3 = simult_sparse_coding(A*D3,A*Y,d,k,1);

%% reconsruct image from distinct patches
Ycomp1 = D1*X1;
YcompImage1 = col2im(Ycomp1, [PATCH.ROWS PATCH.COLS], [size(ima,1) size(ima,2)], 'distinct');

Ycomp2 = D2*X2;
YcompImage2 = col2im(Ycomp2, [PATCH.ROWS PATCH.COLS], [size(ima,1) size(ima,2)], 'distinct');

Ycomp3 = D3*X3;
YcompImage3 = col2im(Ycomp3, [PATCH.ROWS PATCH.COLS], [size(ima,1) size(ima,2)], 'distinct');

% %% display results
figure(1); imshow(uint8(YcompImage1)); title('KSVD');
figure(2); imshow(uint8(YcompImage2));  title('BKSVD-SAC');
figure(3); imshow(uint8(YcompImage3));  title('BKSVD');



%% save results
imwrite(uint8(YcompImage1),['SuperRes-KSVD.jpg']);
imwrite(uint8(YcompImage2),['SuperRes-BKSVD-SAC.jpg']);
imwrite(uint8(YcompImage3),['SuperRes-BKSVD.jpg']);



%%%%%%%%%%%%%%%%%%%%% overlapping patches %%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract all patches in training image
Y = im2col(double(ima), [PATCH.ROWS PATCH.COLS], 'sliding');
[rows,cols] = size(ima);


%% reconstruct signals from compressed vectors
X1 = simult_sparse_coding(A*D1,A*Y,d0,k*s,1);
X2 = simult_sparse_coding(A*D2,A*Y,d2,k,1);
X3 = simult_sparse_coding(A*D3,A*Y,d,k,1);

%% reconsruct patches
Ycomp1 = D1*X1;
Ycomp2 = D2*X2;
Ycomp3 = D3*X3;

%% reconstruct image be averaging over all patches that include each pixel
YcompImage1 = zeros(size(ima));
YcompImage2 = zeros(size(ima));
YcompImage3 = zeros(size(ima));
normalizer = ones(size(ima));
hmin = ceil(PATCH.ROWS/2);
hmax = floor(PATCH.ROWS/2)-1;
wmin = ceil(PATCH.COLS/2);
wmax = floor(PATCH.COLS/2)-1;
g=gausswin(8);
kernel = g*g'; kernel = kernel / sum(sum(kernel));

sum_kernel =  zeros(size(ima));
count = 1;
for j=1+wmin:cols-wmax
    for i=1+hmin:rows-hmax         
        % weighted mean
        YcompImage1(i-hmin:i+hmax,j-wmin:j+wmax) = ...
            YcompImage1(i-hmin:i+hmax,j-wmin:j+wmax) + ...
            reshape(Ycomp1(:,count),PATCH.ROWS,PATCH.COLS).*kernel;
        YcompImage2(i-hmin:i+hmax,j-wmin:j+wmax) = ...
            YcompImage2(i-hmin:i+hmax,j-wmin:j+wmax) + ...
            reshape(Ycomp2(:,count),PATCH.ROWS,PATCH.COLS).*kernel;
        YcompImage3(i-hmin:i+hmax,j-wmin:j+wmax) = ...
            YcompImage3(i-hmin:i+hmax,j-wmin:j+wmax) + ...
            reshape(Ycomp3(:,count),PATCH.ROWS,PATCH.COLS).*kernel;
        
        sum_kernel(i-hmin:i+hmax,j-wmin:j+wmax) = sum_kernel(i-hmin:i+hmax,j-wmin:j+wmax)+kernel;
         count = count+1;
    end
end

%% display results
figure(11); imshow(uint8(YcompImage1./sum_kernel)); title('KSVD with overlap');
figure(12); imshow(uint8(YcompImage2./sum_kernel));  title('BKSVD-SAC with overlap');
figure(13); imshow(uint8(YcompImage3./sum_kernel));  title('BKSVD with overlap');


%% save results
imwrite(uint8(YcompImage1./sum_kernel),['SuperRes-overlap-KSVD.jpg']);
imwrite(uint8(YcompImage2./sum_kernel),['SuperRes-overlap-BKSVD-SAC.jpg']);
imwrite(uint8(YcompImage3./sum_kernel),['SuperRes-overlap-BKSVD.jpg']);







