% In this code the aim is to modify the standard MMCA problem by adding a
% term for learning local dictionaries obtained from patches of the sources
% during the separation process. The following cost function is to
% minimized:
% {alpha_i,S,D,A} = min lambda*||Y-AS||+ lambda_2||s_i*Phi||+
% sum(||si*Rik-Di*alpha_k||+etta*||alpha_k||)
% Phi can be a basis e.g. DCT transform
% The focus here is to deal with noisy cases and also different varying
% reqularization parameters.
% Created  06/05/2011

clc
clear all
close all

nmax= 80; % maximum number of iterations
% add noise
SNR_db = 20;
% load random states for repeatable experiments
% load RandomStates
% rand('state', rand_state);
% randn('state', randn_state);

% im1 = double(imread('lena.png'));
im1 = double(imread('boat256.png'));
% im1 = double(imread('mask_512x512.bmp'));
im2 = double(imread('barbara256.png'));

% im4 = double(imread('paraty256.tif'));
im3 = double(imread('paraty256_2.tif'));
% im3 = double(imread('sosdelrey.tif'));
im4 = double(imread('pakhawaj.tif'));
% im4 = double(imread('texture4.tif'));

% im3 = randn(size(im1))*90;

% % texture
% im1 = double(imread('D1.gif'));
% im2 = double(imread('D2.gif'));
% im3 = double(imread('D3.gif'));
% im4 = double(imread('D4.gif'));
% 

% resize
% im1 = imresize(im1,[128 128]);
% im2 = imresize(im2,[128 128]);
% im3 = imresize(im3,[128 128]);
% im4 = imresize(im4,[128 128]);

NbSources = 4;
NbChannels = 10;
% generate the mixture
A = colnorm(NbChannels,NbSources);
while sum(svd(A)>.2) < min(NbSources,NbChannels)
    A = colnorm(NbChannels,NbSources);
end
% A = orth(A);

X = [reshape(im1,[1 prod(size(im1))]);...
    reshape(im2,[1 prod(size(im1))]);...
                    reshape(im3,[1 prod(size(im1))]);...
                    reshape(im4,[1 prod(size(im1))])...
    ];

Y = A*X;
Y_zeronoise = Y;
Y = Y + std2(Y)*10^(-SNR_db/20)*randn(size(Y)); %add noise type I
sigma = std2(Y_zeronoise)*10^(-SNR_db/20);

% ----------------------  add noise type II --------------------- %
sigmamin = 0.5;
% noise = randn(size(Y)) * sigmamin;
% Y = Y + noise;                

[n] = size(A,1);
[n,N] = size(Y);
[NbChannels,NbSamples] = size(Y);
stepsize = [6 6]; % was 2,2
dicsz = 256;
patchsz = 8;
const=sqrt(1.15);
numIteration = 2;

% initial DCT dictionary, if needed
Dictionary=zeros(patchsz,sqrt(dicsz));
for k=0:1:sqrt(dicsz)-1,
    V=cos([0:1:patchsz-1]*k*pi/sqrt(dicsz));
    if k>0, V=V-mean(V); end;
    Dictionary(:,k+1)=V/norm(V);
end;
Dictionary=kron(Dictionary,Dictionary);
DCT=Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
% save initial DCT dcitioanry for all sources
for ff = 1:NbSources
    dict(:,:,ff) = DCT;
end
Dictionary = DCT;

Y_ = Y;  % store Y
% removing the average
for ll=1:NbChannels
    Y(ll,:) = Y(ll,:) - mean(Y(ll,:))*ones(size(Y(ll,:)));
end

etta = 5; % stepsize for gradient descent in dictionary update
ettastep = (etta-0.1)/nmax;
temp = 0;
his2 = [];his1 =[];svdA = [];his = [];his3=[0;0];his4=[];
indx = 1:dicsz;

% multi-initialization parameteres
params.stepsize = stepsize;
params.patchsz = patchsz;
params.KS = 5;
params.Kmin = 0;
params.orth = 'orth'
params.dicsz = dicsz;
params.NbResets = 5;
params.nmax = 20;
params.const = const;
params.numIteration = 1;
params.lambdaCoef = 5;
params.sigma0=10;
params.NbSources=NbSources;
params.sigmamin=sigmamin;
params.imsize = size(im1);
params.trans = 'dct';

% % %==== 111 multiinitialization to find the best AA
% [AA error_AAoriginalA error_YS SEst_r] = multiInit(Y,NbSources,params,A);
% E = Y-AA*SEst_r;       %residual
% %-----------------------------------------------

%==== 222 Or signle random initilalization
% initial values
AA = colnorm(NbChannels,NbSources); %  randomlise and normalizing the mxing matrix
AA= orth(AA); % build orthonormal basis
SEst_r = AA'*Y;
E = Y-AA*SEst_r;       %residual
%-------------------------------------------------

%================Main parameters for proposed LOCAL MCA====================
% sigma0 = 10;
% sigma = sigma0;
offset = 0;     % offset for sigmamin,
sig = (sigma-(sigmamin+offset))/(nmax);
% lambda=30/sigma^2;
% adjust threshold
lambda = 30/sigma;
% lambda = 5;
% lambda=exp(.5./sigma);
% lambda0 = 2;
% lambda = lambda0;
% lambda_max = 20;              % denoising parameter, minimum allowed value
% lamstep = (lambda_max - lambda)/(nmax);
% mu = muCoef/sigma;              % increasing mu
% lambda = lambda_max*lambda;         % decreasing mu
% mu = 0;
%==========================================================================

%================Main parameters for GLOBAL MMCA===========================
Kmin=0;
KS = 0;
DKS = (KS-Kmin)/(nmax);  %-- Should be changed - Typical values for an appropriate thresholding
% lambda2 = 0;
% lambda2 = sigma0/30-mu;
% qmf = MakeONFilter('Battle',5);
%==========================================================================
% Dict_all = [];

[xn,Q]=nearest_w(AA,A); % solve permutation
% subplot(3,2,5),hold on
% plot(([xn(:,1) A(:,1)]));(([xn(:,1) A(:,1)]))
% subplot(3,2,6),hold on
% plot(([xn(:,2) A(:,2)]));(([xn(:,2) A(:,2)]))


s = 3; %dimension of blocks
k = 4;
K = 256;           %Number of columns in dictionary
max_it = 1;       %Nr of iterations for the algorithm to converge. 
B = K/s;          %Number of blocks in D. Should be an order of magnitude higher than the 
                  %block sparsity k, otherwise the representations wont be sparse
d0 = 1:K;         %block structure with K blocks of size 1 (i.e. block structure ignored)
d = repmat(1:B, s,1); d = d(:)';  %block structure with B blocks of size s
Re_all=[];
atom_count = [];
tic

for (pp=1:nmax)
    
    for ff = 1:NbSources      
        Ei = E+AA(:,ff)*SEst_r(ff,:);
        x = reshape(SEst_r(ff,:),size(im1));
        % Creating image patches
        blocknum = prod(floor((size(x)-[patchsz patchsz])./stepsize) + 1); 
        Data=zeros(patchsz^2,blocknum);
        cnt=1;
        for j=1:stepsize(1):size(x,1)-patchsz+1
            for i=1:stepsize(2):size(x,2)-patchsz+1
                patch=x(i:i+patchsz-1,j:j+patchsz-1);
                Data(:,cnt)=patch(:); % Data has size 64*1, patches is 8*8
                cnt=cnt+1; % 15626
            end;
        end; 
        dc = 0;
        %         [Data, dc] = remove_dc(Data,'columns');
        Dictionary = dict(:,:,ff);
        %------------------- step 2 --- update the dictionary------------
        for jj = 1:numIteration
            %------------------- step 1 --- update the sparse CoefMatrix------------
%             CoefMatrix=omp2(Dictionary'*Data,sum(Data.*Data),Dictionary'*Dictionary,patchsz*sigma*const);
%             All_Coef{ff} = CoefMatrix;
%             if   sum(sum(CoefMatrix~=0)) < 1
%                 error('sigma0 MUST be smaller....');
%             end

%             %=1111======== Gradient Descent =================================
%             Dictionary = Dictionary+etta*(Data-Dictionary*CoefMatrix)*CoefMatrix';
%             Dictionary = colnorm(Dictionary);
%             his4 = [his4 norm(Data-Dictionary*CoefMatrix,'fro')];
%             figure(1);
%             subplot(3,2,6), plot(his4');title('dic learning error');
%             %===============================================================
            %===22222============OR KSVD ========================================
            % Update the dictionary
%             Dictionary(:,1)=Dictionary(:,1); % the DC term remain unchanged
            old_dictionary = Dictionary;
            % --------------------------------------------------------------- %
            [Dictionary,CoefMatrix,d2] = BKSVD(Data,K,k,s,d0,max_it,old_dictionary);
            All_Coef{ff} = CoefMatrix;
%             --------------------------------------------------------------- %
%             for j=2:1:size(Dictionary,2) % implement the K-SVD algorithm
%                 relevantDataIndices=find(CoefMatrix(j,:)); %Find indices of nonzero elements.
%                 if ~isempty(relevantDataIndices)
%                     tmpCoefMatrix=CoefMatrix(:,relevantDataIndices);
%                     tmpCoefMatrix(j,:)=0;
%                     errors=Data(:,relevantDataIndices)-Dictionary*tmpCoefMatrix;
%                     [betterDictionaryElement,singularValue,betaVector]=svds(errors,1);
%                     CoefMatrix(j,relevantDataIndices)=singularValue*betaVector';
%                     Dictionary(:,j)=betterDictionaryElement;
%                 end;
%             end;
            % ---------------------------------------------------- %
        end
        dict(:,:,ff) = Dictionary;
        %================================================================
        
        %------------------- step 4 --- update the sources SEst_r -----------------
        my_y = reshape(AA(:,ff)'*Ei,size(im1));
        
        %         % 111111======> LOCAL ONLY USING LEARNED DCITIONARIES
        %         [y(:,:,ff) yout1]=RecoverImageWithGlobalTrans(my_y,mu,Dictionary,CoefMatrix,blocknum,stepsize,dc);
        %         SEst_r(ff,:) = reshape(y(:,:,ff), [1 NbSamples]);
        %         % ---------------------------------------------------------------
        
        %         %222222=====> GLOBAL MMCA ONLY using soft thresholding
        %         [my_y, dc] = remove_dc(my_y);
        %         temp = FWT2_PO(my_y,1,qmf);
        %         temp = temp(:).*(abs(temp(:)) > KS*mad2(temp(:)));
        %         my_y = IWT2_PO(reshape(temp,size(im1)),1,qmf);
        %         SEst_r(ff,:) = my_y(:);
        %         %---------------------------------------------------------------
        
        %         % 333333======> GLOBAL+LOCAL USING SIGN OPERATION
        %         [y(:,:,ff) yout1]=RecoverImageWithGlobalTransSign(my_y,x,mu,lambda2,Dictionary,CoefMatrix,blocknum,stepsize,dc,'wavelet');
        %         SEst_r(ff,:) = reshape(y(:,:,ff), [1 NbSamples]);
        %         % ---------------------------------------------------------------
        
        %         % 44444 =====> GLOBAL+LOCAL Using Soft Threshold PLUS Local Dictionary
        %         [y(:,:,ff) yout1]=RecoverImageWithGlobalTrans(my_y,mu,Dictionary,CoefMatrix,blocknum,stepsize,dc);
        %         [y(:,:,ff) dc] = remove_dc(y(:,:,ff));
        %         temp = FWT2_PO(y(:,:,ff),1,qmf);
        %         temp = temp(:).*(abs(temp(:)) > KS*mad2(temp(:)));
        %         SEst_r(ff,:) = reshape(IWT2_PO(reshape(temp,size(im1)),1,qmf),[1 prod(size(im1))]);
        %         SEst_r(ff,:) = add_dc(SEst_r(ff,:),dc);
        %         %---------------------------------------------------------------
        %
        
        % 55555 =====> GLOBAL+LOCAL Using Soft Threshold PLUS Local Dictionary
        
        % We use an improved (and more correct) version of the 44444 method
        [y(:,:,ff) yout1 Weight]=RecoverImageWithGlobalTrans3(my_y,lambda,Dictionary,CoefMatrix,blocknum,stepsize,dc);
        %         [y(:,:,ff) dc] = remove_dc(y(:,:,ff));
        
        switch(params.trans)
            case  'wavelet'
                temp = FWT2_PO(y(:,:,ff),1,qmf); % wavelet transform
                temp = temp(:).*(abs(temp(:)) > KS*mad2(temp(:)));
                my_y2 = IWT2_PO(reshape(temp,imsize),1,qmf);
            case 'dct' % curvelet transform
                temp = dct2(y(:,:,ff));
                temp = temp(:).*(abs(temp(:)) > KS*mad2(temp(:)));
                my_y2 = idct2(reshape(temp,size(im1)));
        end
        %         my_y2 = add_dc(my_y2,dc);
        my_y2 = my_y2./(lambda+Weight);
        SEst_r(ff,:) = my_y2(:);
        Ei = E+AA(:,ff)*SEst_r(ff,:);
        % ---------------------------------------------------------------      
%                 figure(7);
%                 if ff == 1
%                     subplot(2,1,1);imshow(my_y2,[]);
%                 else
%                     subplot(2,1,2);imshow(my_y2,[]);
%                 end
%                 pause(.1);
        %------------------- step 3 --- update the mixing matrix columns-----------
        AA(:,ff) = colnorm(Ei*SEst_r(ff,:)');
%         AA2(:,ff)= colnorm(Ei*SEst_r(ff,:)');
        Ei = E+AA(:,ff)*SEst_r(ff,:);
        E = Y-AA*SEst_r;
        his3(ff,pp)=norm(Data-Dictionary*CoefMatrix);
    end
    
    count = nnz(All_Coef{ff});
    atom_count = [atom_count;count];
    %     AA=AA2;    
%     svd(AA)
%     Dict_all = [Dict_all; Dictionary];
    svdA = [svdA svd(AA)]; % returns a vector containing the singular values.
    
    his = [his norm(E,'fro')];
    his1 = [his1 norm(Y-AA*SEst_r,'fro')];
    his2 = [his2 std2(Ei)];
    [xn,Q]=nearest_w(AA,A);
    his4 = [his4 norm(temp(:),1)];
%     if(rem(pp,2)==0)
%         figure(1);
%         subplot(3,2,1);plot([his' his1']);title('||Y-AS||_F');
%         subplot(3,2,2), plot((sum(his3,1)+his)');title('sum of all errors');
%         subplot(3,2,3), plot(sum(his3,1)');title('sum of all errors for DL');
%         subplot(3,2,4),plot([svdA']);title('svd');
%         subplot(3,2,5),hold on
%         plot(([xn(:,1) A(:,1)]));(([xn(:,1) A(:,1)]))
%         subplot(3,2,6),hold on,
%         plot(([xn(:,2) A(:,2)]));(([xn(:,2) A(:,2)]))
%     end
    
%     if sum(svd(AA) < 0.02) > 0
%         AA = orth(AA);
%         E = Y - AA*SEst_r;
%         figure(1);clf;
%     end
    
    %========== GLOBAL MMCA PARAMETERS=======================================
        
%     lambda = lambda + lamstep;
        sigma = sigma - sig;
%         sigma = sigma0*exp(1./nmax*log((sigmamin+offset)/sigma0)*pp);
    %plot(sigma0*exp(1./nmax*log((sigmamin+offset)/sigma0)*(1:nmax)));
     

%     lambda = 30/sigma;
        
    %     lambda=exp(.5./sigma);
    %         mu = 3/sigma;    
    %     lambda = 5;
%         lambda = lambda0*exp(log(lambda_max/lambda0)/nmax*pp);
    %plot(lambda0*exp(log(lambda_max/lambda0)/nmax*(1:nmax)));
        %     etta = etta - ettastep;

    %==========================================================================
    
    
    %========== LOCAL MMCA PARAMETERS=======================================
    %     lambda2 = lambda2- DKS;
    %     lambda2 =  sigma/2;
    %       lambda2 = sigma0/30-mu;
    %     lambda2 = KS;
    KS = KS - DKS;
    %==========================================================================
    E = Y - AA*SEst_r;
    Re = X - SEst_r;
    Re = norm(Re,2)/10;
    Re_all = [Re_all; Re];
    disp(pp)
end
toc

S=SEst_r;
% S = pinv(AA)*Y;

[xn,Q]=nearest_w(AA,A);
% for i=1:NbSources; figure(4);plot(([xn(:,i) A(:,i)]));(([xn(:,i) A(:,i)])),pause,end

S = Q'*S;       % apply permutation

% --------- computing the SNR -----------------
% im_mix1 = mat2gray(reshape(Y_zeronoise(1,:),size(im1)))*255;
% im_mix_noisy1 = mat2gray(reshape(Y(1,:),size(im1)))*255;
% PSNR_mixture = 20*log10(255 * sqrt(numel(im_mix1)) / norm(im_mix1(:)-im_mix_noisy1(:)));
% 
im_source1 = mat2gray(im1)*255;
im_source_noisy1 = mat2gray(reshape(S(1,:),size(im1)))*255;
PSNR_source1 = 20*log10(255 * sqrt(numel(im_source1)) / norm(im_source1(:)-im_source_noisy1(:)));

im_source2 = mat2gray(im2)*255;
im_source_noisy2 = mat2gray(reshape(S(2,:),size(im2)))*255;
PSNR_source2 = 20*log10(255 * sqrt(numel(im_source2)) / norm(im_source2(:)-im_source_noisy2(:)));

im_source3 = mat2gray(im3)*255;
im_source_noisy3 = mat2gray(reshape(S(3,:),size(im3)))*255;
PSNR_source3 = 20*log10(255 * sqrt(numel(im_source3)) / norm(im_source3(:)-im_source_noisy3(:)));

im_source4 = mat2gray(im4)*255;
im_source_noisy4 = mat2gray(reshape(S(4,:),size(im4)))*255;
PSNR_source4 = 20*log10(255 * sqrt(numel(im_source4)) / norm(im_source4(:)-im_source_noisy4(:)));

avg_PNSR_source = (PSNR_source1 + PSNR_source2 + PSNR_source3 + PSNR_source4)/4

R1 = abs(corr2(im1,im_source_noisy1));
R2 = abs(corr2(im2,im_source_noisy2));
R3 = abs(corr2(im3,im_source_noisy3));
R4 = abs(corr2(im4,im_source_noisy4));
avg_R = (R1+R2+R3+R4)/4


% 
% figure(2);
% subplot(2,2,1);imshow(reshape(Y(1,:),size(im1)),[]);title('mixture');
% title(sprintf('PSNR = %.2fdB', PSNR_mixture));
% subplot(2,2,2);imshow(reshape(Y(2,:),size(im1)),[]);title('mixture');
% % subplot(2,2,3);imshow(reshape(Y(3,:),size(im1)),[]);title('mixture');
% % subplot(2,2,4);imshow(reshape(mean(Y,1),size(im1)),[]);title('average mixture');
% 
% figure(3);
% subplot(2,2,1);imshow(reshape(S(1,:),size(im1)),[]);
% title(sprintf('source image, PSNR = %.2fdB', PSNR_source1));
% subplot(2,2,2);imshow(reshape(-S(2,:),size(im1)),[]);title('sources');
% title(sprintf('source image, PSNR = %.2fdB', PSNR_source2));
% subplot(2,2,3);imshow(reshape(S(3,:),size(im1)),[]);title('sources');
% subplot(2,2,4);title('sources');imshow(reshape(S(4,:),size(im1)),[]);
% 
figure;
dictimg = showdict(dict(:,:,1),[1 1]*patchsz,round(sqrt(dicsz)),round(sqrt(dicsz)),'lines','highcontrast');
imshow(imresize(dictimg,2,'nearest'));
title('Trained dictionary');

for imiter = 1:4
    figure;
    imshow(reshape(S(imiter,:),256,256),[]);
end
% mixing matrix criterion
mmc(A,pinv(AA),Q)
figure
semilogy(1:nmax,his1)
% plot(Re_all);

%%

% % Apply ICA to the sparse mixtures
% [S_ica, A_ica, e_s_ica] = fastica(Y,'numOfIC',NbSources);
% 
% figure(6);
% subplot(2,2,1);imshow(reshape(S_ica(1,:),size(im1)),[]);title('sources');
% subplot(2,2,2);imshow(reshape(S_ica(2,:),size(im1)),[]);title('sources');
% subplot(2,2,3);imshow(reshape(S_ica(3,:),size(im1)),[]);title('sources');
% subplot(2,2,4);imshow(reshape(S_ica(4,:),size(im1)),[]);title('sources');
% 
% 
% [xn,Q]=nearest_w(A_ica,A);
% xn=colnorm(xn);
% for i=1:NbSources; figure(4);plot(([xn(:,i) A(:,i)]));(([xn(:,i) A(:,i)])),pause,end
% 
