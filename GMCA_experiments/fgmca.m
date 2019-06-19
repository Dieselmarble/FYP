function [piA,S,his1] = fgmca(X,NbSources,nmax,Kmin)

%  Sped-Up GMCA
%
%  INPUT - X : the sparse data - channels * samples
%		 NbSources : number of sources to estimate
%		 nmax : number of iterations (typically 100)
%		 Kmin : last k-mad threshold (typically 3)
%
%  OUTPUT - piA : the pseudo inverse of the mixing matrix
%		  S   : estimated/thresholded sparse sources
%
%  REQUIRED TOOLBOX :  Wavelab 850 - http://www-stat.stanford.edu/~wavelab/
%  FACULTATIVE TOOLBOXES : MCALab - http://www.greyc.ensicaen.fr/~jfadili/demos/WaveRestore/downloads/mcalab.html
%					CurveLab - http://www.curvelet.org/
%
%  Version : 7th of August 2007
%		   CEA-Saclay / DAPNIA/SEDI-Sap
%		   J.Bobin
%
%  

[NbChannels,NbSamples] = size(X);
his1 = [];
for ll=1:NbChannels
	
	X(ll,:) = X(ll,:) - mean(X(ll,:))*ones(size(X(ll,:)));

end 

KS = 15;
DKS = (KS-Kmin)/nmax;  %-- Should be changed - Typical values for an appropriate thresholding

AA = randn(NbChannels,NbSources);

SEst_r = zeros(NbSources,NbSamples);

for pp=1:nmax
            
            piA = inv(AA'*AA)*AA'; % definitaion of pesudo inverse
            
            SEst_r = piA*X; % estimation of sources
            
            SigmaSources = zeros(1,NbSources); 
            
            for ff = 1:NbSources        
                   SEst_r(ff,:) = SEst_r(ff,:).*(abs(SEst_r(ff,:)) > KS*mad(SEst_r(ff,:))); % soft thresholding operater
%                    SEst_r(ff,:) = sign(SEst_r(ff,:)).*max(abs(SEst_r(ff,:)) - KS/2,0);%hard thresholding
                   SigmaSources(ff) = std(SEst_r(ff,:)); % choose the support sets

            end
            
            indd = find(SigmaSources > 1e-9);
            
            if length(indd) > 0
            
                    AA(:,indd) = X*SEst_r(indd,:)'*inv(SEst_r(indd,:)*SEst_r(indd,:)'); %least quare estimation of mixing matrix A
                    
                    for ff = 1:length(indd)
                    
                    	AA(:,indd(ff)) = AA(:,indd(ff))/norm(AA(:,indd(ff))); %normalise matrix A
                                                       
                    end
             end
                           
             KS = KS - DKS; % decrease threshold 
             his1 = [his1 norm(X-AA*SEst_r,'fro')];
end

S = SEst_r;
