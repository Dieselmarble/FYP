% nGMCAhard_alg.m - This file is part of nGMCALab.
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
% Aims at solving argmin_(A >= 0, S >= 0) 1 / 2 * ||Y - A * S||_2^2 + lambda * ||S||_0
% using iterative hard-thresholding.
% Parameter lambda is decreasing during the iterations and set to
% tau_MAD*sigma_noise at the end of the algorithm, where tau_MAD is a
% constant (preferably in [1,3]) and sigma_noise is an online estimate of
% the noise standard deviation.
%
% For more information, this algorithm is described as nGMCA^H in:
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
%   and refinement phase in pourcent of the iterations (default: 0.80)
% - A.MaximumIteration and S.MaximumIteration: maximum number of iteration
%   of the Forward-Backward subroutines in A and S (default: 80).
% - A.RelativeDifferenceTolerance and A.RelativeDifferenceTolerance:
%   relative difference tolerance for convergence of the Forward-Backward
%   sub-iterations in A and S (default: 0.00001).
% - S.tau_MAD: constant coefficient for the final threshold
%   computation (default: 1). 
% - lambdaInf: in order to predefine the final lambda (not advised,
%   except for underdetermined settings).
%
% %Output: structure with fields:
% - lambda: the mean final lambda.
% - A: the mixing matrix.
% - S: the sparse sources matrix.
% With the approximation Y \approx A*S.
%
% %Example of use:
% m = 200; %number of observations
% n = 200; %number of source samples
% r = 5; %number of sources
% ref_A = rand(m,r); %non-negative mixtures
% ref_S = rand(r,m) .* (rand(r, m) > 0.9); %non-negative and sparse sources
% Y = ref_A * ref_S + 0.1*randn(m, n); %simulated data
%
% parameters.Y = Y;
% parameters.rank = r;
% algo = nGMCAhard_alg();
% result = ApplyAlgorithmn(algo, parameters); %perform decomposition
% 
% plot(result.S');
%
function algorithm = nGMCAhard_alg(parameters)
    if nargin < 1
        parameters = [];
    end
    algorithm.parameters = parameters;

    algorithm.initialize = @m_algo_initialize;
    algorithm.iterate = @m_algo_iterate;
    algorithm.terminate = @m_algo_terminate;

    algorithm.name = 'nGMCA^H';
end



function [data,parameters] = m_algo_initialize(data,parameters)
% initialization of the data structure (and optionnally the function value)

%% check for required parameters
if ~isfield(parameters, 'Y')
    error('In nGMCA_alg, parameter field Y (input data) is required.\n')
end
if ~isfield(parameters, 'rank')
    error('In nGMCA_alg, parameter field rank (number of sources) is required.\n')
end

name = 'nGMCAhard_alg';
parameters = addMissingParameter(parameters, 'verbose', 0, name);
parameters = addMissingParameter(parameters, 'MaximumIteration', 500, name);
parameters = addMissingParameter(parameters, 'phaseRatio', 0.80, name);

parameters = addMissingParameter(parameters, 'S.tau_MAD', 1, name);
parameters = addMissingParameter(parameters, 'S.MaximumIteration', 80, name);
parameters = addMissingParameter(parameters, 'S.RelativeDifferenceTolerance', 0.00001, name);
parameters = addMissingParameter(parameters, 'A.MaximumIteration', 80, name);
parameters = addMissingParameter(parameters, 'A.RelativeDifferenceTolerance', 0.00001, name);

%% random initialization if not given
if ~isfield(data, 'A')
    data = NMFwarmInit(parameters.Y, parameters.rank);
end

%% prevent bad initialization (happening for low number of observations)%%%%%%%%%%%%%%%%%%%%%%%%%
data = NMFnormalization(data, 'A');
data = reinitializeNullSources(data, parameters.Y, parameters.verbose);


%% first threshold
data = NMFnormalization(data, 'A');
data.doNotCheckConvergence = 1;
parameters.refinement_beginning = floor(parameters.phaseRatio * parameters.MaximumIteration);

end






function [data, relativeDifference] = m_algo_iterate(data, parameters)
% iterate and update the data (and optionnally the function value)

    %prepare variables
    prevAS = [al(data.A); al(data.S)];
    n = size(data.S, 2);
    
    %% update thresholds
    data = m_thresholdManagement(data, parameters);

    %% update S
    data.S = nonnegativeSparseUpdate_hard_S(data.A, data.S, parameters.Y,...
        (data.kmad .* data.sigs) * ones(1,n), parameters.S);

    %reinitialize lines of S if need be
    data = reinitializeNullSources(data, parameters.Y, parameters.verbose);
        
    
    %% update A
    %preconditioning
    data = NMFnormalization(data, 'S');
    %Constraining the columns norm does not change anything in this case since the
    %problem is not convex.

    %update
	data.A = nonnegativeSparseUpdate(data.A', data.S * data.S',...
        data.S * parameters.Y', 0, parameters.A)';

    
    %% update lambda
    data = NMFnormalization(data, 'A');
    
    
    %% relative difference
    relativeDifference = norm([al(data.A); al(data.S)] - prevAS) / norm(prevAS);

end


function result = m_algo_terminate(data)
% output the result of the algorithm

result.A = data.A;
result.S = data.S;
result.lambda = mean(data.kmad .* data.sigs);

end


%% ALGORITHM SPECIFIC FUNCTION
%%


function S = nonnegativeSparseUpdate_hard_S(A, S0, Y, lambda, parameters)
%aims at solving (non-convex) problem
%argmin_(S >= 0) 1 / 2*||Y - A * S||_2^2 + lambda*|| S ||_0
%with starting point S0

%precomputation
H = A' * A;
AtY = A' * Y;
L = max(eig(H));

%initializations
S = S0;
prev_S = S;
t = 1;

%operators
gradient = @(s) H * s - AtY;
proximal = @(x, threshold) x .* (x > threshold);
al = @(x) x(:);

for k = 1 : parameters.MaximumIteration
    t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
    w = (t - 1) / t_next;
    t = t_next;
    
    R = (1 + w) * S - w * prev_S;
    prev_S = S;
    S = proximal(R - gradient(R) / L, lambda / L);

    if norm(al(prev_S - S), 'fro') / norm(al(S), 'fro') < parameters.RelativeDifferenceTolerance
        break
    end
end

end


%% Parameter update
function data = m_thresholdManagement(data, parameters)
%the strategy with hard-thresholding must differ from the one with
%soft-thresholding: in order to separate correctly, the threshold must be
%at the level of the sources in the beginning. The threshold ultimately 
%applying on the gradient, it must be of the level of the gradient in the
%end.

    %precomputations
    n = size(data.S, 2);
    H = data.A' * data.A;
    tau = 1 / max(eig(H));
    M = floor(parameters.phaseRatio * parameters.MaximumIteration - data.iteration);
    ratio = data.iteration / (parameters.phaseRatio * parameters.MaximumIteration);%ratio of activation
    
    %S with one step of gradient (divided by tau)
    updS = abs(data.S / tau -(H * data.S - data.A' * parameters.Y));
    
    %estimation of sigmas on the signals
    sigs = dimMADstd(updS, 2);
    %estimation of sigma on the residual
    Rsig = dimMADstd(al(parameters.Y - data.A * data.S), 1);
    %weighting
    data.sigs = Rsig * min(1, ratio) + max(1 - ratio, 0) * sigs;

    if M>0
        updS = al(updS ./ (data.sigs * ones(1, n)));
        updS = sort(updS, 1, 'descend');
        N = sum(updS > parameters.S.tau_MAD);
        data.kmad = max(parameters.S.tau_MAD, updS(min(N, max(1, 1 + floor(ratio * N)))));
    else
        data.kmad = parameters.S.tau_MAD;
    end

end
