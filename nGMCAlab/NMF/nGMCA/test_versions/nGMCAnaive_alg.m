% nGMCAnaive_alg.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13, last modified on 16/7/2014
%
% This software is governed by the CeCILL  license under French law and
% abiding by the rules of distribution of free software.  You can  use,
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info".
%
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability.
%
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or
% data to be ensured and,  more generally, to use and operate it in the
% same conditions as regards security.
%
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.
%
%
% non-negative Generalized Morphological Component Analysis
%
% Aims at solving argmin_(A>=0,S>=0) 1/2*||Y-AS||_2^2+lambda*||S||_0
% using a naive least square and hard-thresholding approach.
% Parameter lambda is decreasing during the iterations and set to
% tau_MAD*sigma_noise at the end of the algorithm, where tau_MAD is a
% constant (preferably in [1,3]) and sigma_noise is an online estimate of
% the noise standard deviation in the sources.
%
% For more information, this algorithm is described as nGMCA^naive in:
% J. Rapin, J.Bobin, A. Larue and J.L. Starck,
% Sparse and Non-negative BSS for Noisy Data,
% IEEE Transactions on Signal Processing, 2013.
% Please use the above reference if using this code in a publication.
%
%
% This function creates a structure triggered by the function ApplyAlgorithm.
% Inputs can be provided through a parameter structure at the creation of the
% algorithm (calling this function) or when using it (calling ApplyAlgorithm).
%
% Required fields of the parameter structure (when calling ApplyAlgorithm if
% not set before):
% - Y: data matrix.
% - rank: number of sources.
%
% Optional fields:
% - verbose: display information (default:1) or not (0).
% - MaximumIteration: number of iterations of the main loop
%   (default: 500, should be enough for convergence).
% - phaseRatio: transition between decreasing thresholding phase
%   and refinement phase in pourcent of the iterations (default: 0.80).
% - tau_MAD: constant coefficient for the final threshold
%   computation (default: 1).
%
% %Output: structure with fields:
% - threshold: the mean final threshold.
% - A: the mixing matrix.
% - S: the sparse sources matrix.
% With the approximation Y \approx A*S.
%
% %Example of use:
% m = 200; %number of observations
% n = 200; %number of source samples
% r = 5; %number of sources
% ref_A=rand(m,r); %non-negative mixtures
% ref_S=rand(r,m).*(rand(r,m)>0.9); %non-negative and sparse sources
% Y = ref_A*ref_S + 0.1*randn(m,n); %simulated data
%
% parameters.Y=Y;
% parameters.rank=r;
% algo=nGMCAnaive_alg();
% result = ApplyAlgorithmn(algo,parameters); %perform decomposition
%
% plot(result.S');
%
function algorithm = nGMCAnaive_alg(parameters)
if nargin<1
    parameters=[];
end
algorithm.parameters=parameters;

algorithm.initialize=@m_algo_initialize;
algorithm.iterate=@m_algo_iterate;
algorithm.terminate=@m_algo_terminate;

algorithm.name='nGMCA^{naive}';
end



function [data,parameters]=m_algo_initialize(data,parameters)
% initialization of the data structure (and optionnally the function value)

%% check for required parameters
if ~isfield(parameters,'Y')
    error('In nGMCA_alg, parameter field Y (input data) is required.\n')
end
if ~isfield(parameters,'rank')
    error('In nGMCA_alg, parameter field rank (number of sources) is required.\n')
end

name = 'nGMCAnaive_alg';
parameters = addMissingParameter(parameters, 'verbose', 0, name);
parameters = addMissingParameter(parameters, 'MaximumIteration', 500, name);
parameters = addMissingParameter(parameters, 'phaseRatio', 0.98, name);%refinement is not helpful in this case
parameters = addMissingParameter(parameters, 'S.tau_MAD', 1, name);


%% random initialization if not given
if ~isfield(data,'A')
    data=NMFwarmInit(parameters.Y,parameters.rank);
end

%% prevent bad initialization (happening for low number of observations)%%%%%%%%%%%%%%%%%%%%%%%%%
data = NMFnormalization(data,'A');
data = reinitializeNullSources(data,parameters.Y,parameters.verbose);


%% first threshold
data=NMFnormalization(data,'A');
data.doNotCheckConvergence=1;
parameters.refinement_beginning=floor(parameters.phaseRatio*parameters.MaximumIteration);

end






function [data,relativeDifference]=m_algo_iterate(data,parameters)
% iterate and update the data (and optionnally the function value)

%prepare variables
prevAS=[al(data.A);al(data.S)];
M=parameters.refinement_beginning-data.iteration;
n=size(data.S,2);

%% computation of the treshold
ind=(sum(data.A)>0);
AdY=data.A(:,ind)\parameters.Y;
ind2=~isnan(sum(AdY,2)); %avoid bad conditioning issues
if sum(ind2)~=length(ind2)
    AdY=AdY(ind2,:); %use only correct sources
    ind=find(ind);
    ind=ind(ind2);
end 


if ~isempty(AdY)
    data.sigs=1.4826*median(abs(AdY),2); %estimate of the SOURCES noise std
    if M>0
        kAdY=al(AdY./(data.sigs*ones(1,n)));
        y=sort(kAdY,1,'descend');
        ratio=data.iteration/parameters.refinement_beginning;
        N=sum(y>parameters.S.tau_MAD); %number of samples above the minimum threshold
        data.kmad=max(parameters.S.tau_MAD,y(min(N,max(1,1+floor(ratio*N)))));
        
    else
        data.kmad=parameters.S.tau_MAD;
        data.doNotCheckConvergence=0;
    end

    % update
    data.S(ind,:)=nonnegativeHardThresholding(AdY,data.kmad*data.sigs*ones(1,n));
end


%% reinitialize lines of S if need be
data = reinitializeNullSources(data,parameters.Y,parameters.verbose);


%% update A
data = NMFnormalization(data,'S');%preconditioning
ind=(sum(data.S,2)>0);
data.A(:,ind)=nonnegativeHardThresholding(parameters.Y/data.S(ind,:),0);



%% relative difference
data=NMFnormalization(data,'A');
relativeDifference=norm([al(data.A);al(data.S)]-prevAS)/norm(prevAS);

end


function result=m_algo_terminate(data)
% output the result of the algorithm

result.A=data.A;
result.S=data.S;
result.threshold=mean(data.kmad.*data.sigs);

end





