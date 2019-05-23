% synthesis_conv_nGMCA_alg.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
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
% non-negative Generalized Morphological Component Analysis with synthesis
% sparse regularization in a non-negative convolutive domain
%
% Aims at solving :
% argmin_(A >= 0, S >= 0) 1 / 2 * || Y - A * S_w * W||_2^2 + lambda * ||S_w||_1
% using proximal methods.
% Parameter lambda is decreasing during the iterations and set to
% tau_MAD * sigma_noise at the end of the algorithm, where tau_MAD is a
% constant (preferably in [1,3]) and sigma_noise is an online estimate of
% the noise standard deviation.
%
% For more information, this algorithm is described as analysis nGMCA in:
% J. Rapin, J.Bobin, A. Larue and J.L. Starck,
% NMF with Sparse Regularizations in Transformed Domains, SIAM Journal on 
% Imaging Sciences, 2014
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
% Optional subfields for S:
% - S.MaximumIteration: number of iterations of the update of S  (default: 24)
% - S.RelativeDifferenceTolerance: relative difference tolerance
%   for convergence of the Forward-Backward sub-iterations (default: 10^-5)
% - S.tau_MAD: constant coefficient for the final threshold computation (default: 2).
% Optional subfields for A:
% - S.filter: non-negative filter to be used for the convolutive transform
% (default: []).
% - S.LaplacianFilterWidth:  if no filter is specified, width of the Laplacian
% filter which will be used (default: 4).
% Optional subfields for A:
% - A.MaximumIteration: number of iterations of the update of S  (default: 80)
% - A.RelativeDifferenceTolerance: relative difference tolerance
%   for convergence of the Forward-Backward sub-iterations (default: 10^-5)

%
% %Output: structure with fields:
% - lambda: the final lambda.
% - A: the mixing matrix.
% - S: the sparse sources matrix.
% With the approximation Y \approx A*S.

function algorithm = synthesis_conv_nGMCA_alg(parameters)

if nargin < 1
    parameters = [];
end

algorithm.parameters = parameters;

algorithm.initialize = @m_algo_initialize;
algorithm.iterate = @m_algo_iterate;
algorithm.terminate = @m_algo_terminate;

algorithm.name = 'syn. nGMCA (conv.)';

end

%% ----------------------------------------------------------------------%%


function [data, parameters] = m_algo_initialize(data, parameters)
% initialization of the data structure

%% check for required parameters
if ~isfield(parameters,'Y')
    error('You need to input the matrix Y to factorize in the parameter structure.\n')
end
if ~isfield(parameters,'rank')
    error('You need to input the number n of signals to identify in the parameter structure.\n')
end


%% optional fields
name = 'synthesis_conv_nGMCA_alg';
%main loop
parameters = addMissingParameter(parameters, 'verbose', 0, name);
verbose = parameters.verbose;
parameters = addMissingParameter(parameters, 'phaseRatio',0.85, name, verbose);
parameters = addMissingParameter(parameters, 'MaximumIteration', parameters.rank * 10, name, verbose);

%A update
parameters = addMissingParameter(parameters, 'A.MaximumIteration', 80, name, verbose);
parameters = addMissingParameter(parameters, 'A.RelativeDifferenceTolerance', 0.00001, name, verbose);
%S update
parameters = addMissingParameter(parameters, 'S.MaximumIteration', 80, name, verbose);
parameters = addMissingParameter(parameters, 'S.filter', [], name, 0);
parameters = addMissingParameter(parameters, 'S.LaplacianFilterWidth', 4, name, 0);
parameters = addMissingParameter(parameters, 'S.RelativeDifferenceTolerance', 0.00001, name, verbose);
parameters = addMissingParameter(parameters, 'S.tau_MAD', 3, name, verbose);


%% random initialization if not provided
if ~isfield(data, 'A')
    data = NMFwarmInit(parameters.Y, parameters.rank);
end


%% filter
if isempty(parameters.S.filter)

    if isfield(parameters,'display')
        fprintf(1, 'Initializing Laplacian filter with width %i.\n', parameters.S.LaplacianFilterWidth);
    end
    
    N = size(parameters.Y, 2);
    half_N = floor(N / 2);
    t = 1 : N;
    f = exp(-abs(t - half_N) / (parameters.S.LaplacianFilterWidth * 0.5 / log(2)));
    f = f / norm(f);
    f = f([half_N : N, 1 : half_N - 1]);
    parameters.S.filter = f;
else
    if sum(parameters.S.filter < 0) > 0
       error('All coefficients of the filter should be non-negative');
    end
    parameters.S.filter = parameters.S.filter / norm(parameters.S.filter);
    ind = find(parameters.S.filter == max(parameters.S.filter));
    ind = ind(1);
    parameters.S.filter = parameters.S.filter([ind : end, 1 : ind - 1]);
end

[gradient, lipschitzConstant, convolution, convolution_T] =...
    convolutiveGradient1D(data.A' * data.A, data.A' * parameters.Y, parameters.S.filter);

data.param_S.gradient = gradient;
data.param_S.LipschitzConstant = lipschitzConstant;
data.param_S.convolution = convolution;
data.param_S.convolution_T = convolution_T;
data.param_S.proximal = @nonnegativeSoftThresholding;
data.param_S.MaximumIteration = parameters.S.MaximumIteration;
data.param_S.RelativeDifferenceTolerance = parameters.S.RelativeDifferenceTolerance;


%% initizalize data.Sw
data.param_S.MaximumIteration= parameters.S.MaximumIteration;
data.param_S.initialization.x = 0 * data.S;
data.param_S.lambda = 0;
result = ApplyAlgorithm(ForwardBackward_alg(),data.param_S);
data.Sw = result.x;


%% first threshold
data = NMFnormalization(data, 'A');
g = data.param_S.convolution_T(data.A' * (parameters.Y - data.A * data.S));
data.lambda = max(abs(g(:))) * ones(size(data.S, 1), 1);

%% other variables
parameters.refinement_beginning = floor(parameters.phaseRatio * parameters.MaximumIteration);
data.doNotCheckConvergence = 1;

end


%% ----------------------------------------------------------------------%%

function [data, relativeDifference] = m_algo_iterate(data, parameters)
% iterate and update the data
prevAS = [al(data.A); al(data.S)];

%% update S
data.param_S.initialization.x = data.Sw;
data.param_S.lambda = data.lambda * ones(1, size(data.S, 2));
result = ApplyAlgorithm(ForwardBackward_alg(), data.param_S);
data.Sw = result.x;
data.S = max(0, data.param_S.convolution(data.Sw));

%% reinitialize lines of S if need be
data = reinitializeNullSources(data, parameters.Y, parameters.verbose);


%% update A
%A column norms must be non-increasing in order to converge to a
%stationary point
%this is only helpful in the refinement steps (since the cost function
%is then fixed)
if data.iteration >= parameters.refinement_beginning
    parameters.A.normConstrained = 1;
    data.doNotCheckConvergence = 0; %in refinement phase, convergence can be checked
end
data = NMFnormalization(data, 'S');


SSt = data.S * data.S';
SYt = data.S * parameters.Y';
data.A = nonnegativeSparseUpdate(data.A', SSt,...
     SYt, 0, parameters.A)';
data = NMFnormalization(data, 'A');

%% update gradient S
[data.param_S.gradient, data.param_S.LipschitzConstant] =...
    convolutiveGradient1D(data.A' * data.A, data.A' * parameters.Y, parameters.S.filter);

%% relative difference
relativeDifference = norm([al(data.A); al(data.S)] - prevAS)/norm(prevAS);

%% Update thresholds
data = m_thresholdManagement(data, parameters);

end






%% ----------------------------------------------------------------------%%

function result = m_algo_terminate(data)
% output the result of the algorithm

result.lambda = data.lambda;
result.A = data.A;
result.S = data.S;
end







%% ALGORITHM SPECIFIC FUNCTION
%%

%% Parameter update
function data = m_thresholdManagement(data,parameters)


if isfield(parameters,'sigmaNoise')
    sigma_grad = parameters.sigmaNoise;
else
    res = data.param_S.convolution_T(parameters.Y - data.A * data.S);
    data.sigma_res = dimMADstd(res, 2);
    sigma_grad = sqrt(sum(bsxfun(@times, data.sigma_res.^2, (data.A.^2)))');
end
data.sig = max(sigma_grad);
if data.iteration <= parameters.refinement_beginning
    data.lambda = min(data.lambda, max(parameters.S.tau_MAD * sigma_grad, data.lambda...
        - 1 / (parameters.refinement_beginning + 1 - data.iteration)...
        * (data.lambda - parameters.S.tau_MAD * sigma_grad)));
    data.doNotCheckConvergence = 1;
else
    data.doNotCheckConvergence = 0;
end

end




%%
%%



