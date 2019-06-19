function [Aout error_AAoriginalA error_YS SEst_r] = multiInit(Y,NbSources,params,originalA)
% In this code the aim is to do a set of multi-initilization using
% different methods to find a good initila value for the main separation
% task
% Created 03/05/2011

[NbChannels NbSamples] = size(Y);
nmax=params.nmax;
NbResets = params.NbResets;
stepsize = params.stepsize;
dicsz = params.dicsz;
patchsz = params.patchsz;
const=params.const;
numIteration = params.numIteration; % number of iteration for DL part
imsize = params.imsize;
trans = params.trans;


Y_ = Y;  % store Y
% % removing the average
% for ll=1:NbChannels
%     Y(ll,:) = Y(ll,:) - mean(Y(ll,:))*ones(size(Y(ll,:)));
% end
%

for qq = 1:NbResets
    
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
    
    % initial values
    AA = colnorm(NbChannels,NbSources);
    if (isfield(params,'orth'))
        AA= orth(AA);
    end
    
%     SEst_r = mat2gray(randn(NbSources, NbSamples))*255;
        SEst_r = AA'*Y;
    
    
    E = Y-AA*SEst_r;       %residual
    
    his2 = [];his1 =[];svdA = [];his = [];his3=[0;0];his4=[];
    
    %================Main parameters for proposed LOCAL MCA====================
    sigma0 = params.sigma0;
    sigmamin = params.sigmamin;
    sigma = sigma0;
    sig = (sigma-sigmamin)/(nmax);
    % lambda=30/sigma^2;
    % lambda = 10/(sigma);
    % lambda = 5;
    % lambda=exp(.5./sigma);
    mu = params.muCoef*sigma/30;
    % mu = 0;
    %lambdamax = 30;
    %==========================================================================
    
    
    %================Main parameters for GLOBAL MMCA===========================
    Kmin=params.Kmin;
    KS = params.KS;
    DKS = (KS-Kmin)/(nmax);  %-- Should be changed - Typical values for an appropriate thresholding
    % lambda2 = 0;
    % lambda2 = sigma0/30-mu;
    qmf = MakeONFilter('Battle',5);
    %==========================================================================
    blocknum = prod(floor((imsize-[patchsz patchsz])./stepsize) + 1);
    CoefMatrix=sparse(zeros(dicsz,blocknum));
    clf;
    for (pp=1:nmax)
        
        for ff = 1:NbSources
            
            Ei = E+AA(:,ff)*SEst_r(ff,:);
            x = reshape(SEst_r(ff,:),imsize);
            
            blocknum = prod(floor((size(x)-[patchsz patchsz])./stepsize) + 1);
            Data=zeros(patchsz^2,blocknum);
            cnt=1;
            for j=1:stepsize(1):size(x,1)-patchsz+1
                for i=1:stepsize(2):size(x,1)-patchsz+1
                    patch=x(i:i+patchsz-1,j:j+patchsz-1);
                    Data(:,cnt)=patch(:);
                    cnt=cnt+1;
                end;
            end;
            
            dc = 0;
            %             [Data, dc] = remove_dc(Data,'columns');
            Dictionary = dict(:,:,ff);
            
            for jj = 1:numIteration
                %------------------- step 1 --- update the sparse CoefMatrix------------
                CoefMatrix=omp2(Dictionary'*Data,sum(Data.*Data),Dictionary'*Dictionary,patchsz*sigma*const);
                %                 if   sum(sum(CoefMatrix~=0)) < 1
                %                     error('sigma0 MUST be smaller....');
                %                 end
                %
                %------------------- step 2 --- update the dictionary ---------------------
                %         %========= Gradient Descent =================================
                %         Dictionary = Dictionary+etta*(Data-Dictionary*CoefMatrix)*CoefMatrix';
                %         Dictionary = colnorm(Dictionary);
                %         his4 = [his4 norm(Data-Dictionary*CoefMatrix,'fro')];
                %         figure(1);
                %         subplot(3,2,6), plot(his4');title('dic learning error');
                %         %===============================================================
                %===============OR KSVD ========================================
                % Update the dictionary
                Dictionary(:,1)=Dictionary(:,1); % the DC term remain unchanged
                for j=2:1:size(Dictionary,2)
                    relevantDataIndices=find(CoefMatrix(j,:));
                    if ~isempty(relevantDataIndices)
                        tmpCoefMatrix=CoefMatrix(:,relevantDataIndices);
                        tmpCoefMatrix(j,:)=0;
                        errors=Data(:,relevantDataIndices)-Dictionary*tmpCoefMatrix;
                        [betterDictionaryElement,singularValue,betaVector]=svds(errors,1);
                        CoefMatrix(j,relevantDataIndices)=singularValue*betaVector';
                        Dictionary(:,j)=betterDictionaryElement;
                    end;
                end;
            end
            dict(:,:,ff) = Dictionary;
            %================================================================
            
            %------------------- step 4 --- update the sources SEst_r -----------------
            my_y = reshape(AA(:,ff)'*Ei,imsize);
            
            %         % 111111======> LOCAL ONLY USING LEARNED DCITIONARIES
            %         [y(:,:,ff) yout1]=RecoverImageWithGlobalTrans(my_y,mu,Dictionary,CoefMatrix,blocknum,stepsize,dc);
            %         SEst_r(ff,:) = reshape(y(:,:,ff), [1 NbSamples]);
            %         % ---------------------------------------------------------------
            
            %         %222222=====> GLOBAL MMCA ONLY using soft thresholding
            %         [my_y, dc] = remove_dc(my_y);
            %         temp = FWT2_PO(my_y,1,qmf);
            %         temp = temp(:).*(abs(temp(:)) > KS*mad2(temp(:)));
            %         my_y = IWT2_PO(reshape(temp,imsize),1,qmf);
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
            %         SEst_r(ff,:) = reshape(IWT2_PO(reshape(temp,imsize),1,qmf),[1 prod(imsize)]);
            %         SEst_r(ff,:) = add_dc(SEst_r(ff,:),dc);
            % ---------------------------------------------------------------
            
            
            % 55555 =====> GLOBAL+LOCAL Using Soft Threshold PLUS Local Dictionary
            % We use an improved (and more correct) version of the 44444 method
            [y(:,:,ff) yout1 Weight]=RecoverImageWithGlobalTrans2(my_y,mu,Dictionary,CoefMatrix,blocknum,stepsize,dc);
            %         [y(:,:,ff) dc] = remove_dc(y(:,:,ff));
            
            switch(params.trans)
                case  'wavelet'
                    temp = FWT2_PO(y(:,:,ff),1,qmf);
                    temp = temp(:).*(abs(temp(:)) > KS*mad2(temp(:)));
                    my_y2 = IWT2_PO(reshape(temp,imsize),1,qmf);
                case 'dct'
                    temp = dct2(y(:,:,ff));
                    temp = temp(:).*(abs(temp(:)) > KS*mad2(temp(:)));
                    my_y2 = idct2(reshape(temp,imsize));
            end
            %         my_y2 = add_dc(my_y2,dc);
            my_y2 = my_y2./(1+mu*Weight);
            SEst_r(ff,:) = my_y2(:);
            Ei = E+AA(:,ff)*SEst_r(ff,:);
            % ---------------------------------------------------------------
            
            %         figure(7);subplot(2,1,1);imshow(y(:,:,1),[]);
            %         if pp>1
            %             subplot(2,1,2);imshow(y(:,:,2),[]);pause;
            %         end
            
            %------------------- step 3 --- update the mixing matrix columns-----------
            AA(:,ff) = colnorm(Ei*SEst_r(ff,:)');
            %         AA2(:,ff)= colnorm(Ei*SEst_r(ff,:)');
            Ei = E+AA(:,ff)*SEst_r(ff,:);
            E = Y-AA*SEst_r;
            his3(ff,pp)=norm(Data-Dictionary*CoefMatrix);
        end
        
        svd(AA)
        svdA = [svdA svd(AA)];
        
        his = [his norm(E,'fro')];
        his1 = [his1 norm(Y-AA*SEst_r,'fro')];
        his2 = [his2 std2(Ei)];
        [xn,Q]=nearest_w(AA,originalA);
        his4 = [his4 norm(temp(:),1)];
        
        if(rem(pp,2)==0)
            figure(1);
            subplot(3,2,1);plot([his' his1']);
            subplot(3,2,3), plot(sum(his3,1)');title('sum of all errors for DL');
            subplot(3,2,4),plot([svdA']);title('svd');
            subplot(3,2,5),hold on
            plot(([xn(:,1) originalA(:,1)]));(([xn(:,1) originalA(:,1)]))
            subplot(3,2,6),hold on,
            plot(([xn(:,2) originalA(:,2)]));(([xn(:,2) originalA(:,2)]))
        end
        
        if sum(svd(AA) < 0.02) > 0
            AA = orth(AA);
            E = Y - AA*SEst_r;
            figure(1);clf;
        end
        
        %========== GLOBAL MMCA PARAMETERS=======================================
%         sigma = sigma - sig;
        %     sigma = sigma0*exp(1./nmax*log(sigmamin/sigma0)*pp);
        %     etta = etta - ettastep;
        %     lambda = 30/sigma^2;
        %     lambda=exp(.5./sigma);
%         mu = params.muCoef*sigma/30;
        %     lambda = 5;
        %     lambda = lambda0*exp(log(lambdamax/lambda0)/nmax*pp);
        %==========================================================================
        
        
        %========== LOCAL MMCA PARAMETERS=======================================
        %     lambda2 = lambda2- DKS;
        %     lambda2 =  sigma/2;
        %       lambda2 = sigma0/30-mu;
        %     lambda2 = KS;
        KS = KS - DKS;
        %==========================================================================
        
        disp(pp)
    end
    
    bestA_set(:,:,qq) = AA;
    %     % % removing the average
    %     for ll=1:NbChannels
    %         Y2(ll,:) = Y(ll,:) - mean(Y(ll,:))*ones(size(Y(ll,:)));
    %     end
    error_YS(qq) = norm(Y-AA*SEst_r,'fro');
    error_AAoriginalA(qq) = norm(xn-originalA);
    
end

indx1 = find(error_AAoriginalA == min(error_AAoriginalA));
indx2 = find(error_YS == min(error_YS));

if indx1 == indx2
    disp('!!!!!!! TO INDICES ARE EQUAL !!!!!!!!!!!!!!!');
end
Aout = bestA_set(:,:,indx2);
close all
end