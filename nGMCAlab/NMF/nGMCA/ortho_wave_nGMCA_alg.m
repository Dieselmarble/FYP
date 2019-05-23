% ortho_wave_nGMCA_alg.m - This file is part of nGMCALab.
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
% non-negative Generalized Morphological Component Analysis with sparse 
% regularization in an orthonormal wavelet domain
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
% - S.tau_MAD: constant coefficient for the final threshold computation (default: 2).
% - S.directSparsity: penalize the coarse scale or not (default: 1)
% - S.qmf: wavelet coefficients (use MakeONFilter from the Wavelab toolbox,
% default: Symmlet, 4)
% - S.L: number of scales of the wavelets (default: 3)
% Optional subfields for A:
% - A.MaximumIteration: number of iterations of the update of S  (default: 80)
% - A.RelativeDifferenceTolerance: relative difference tolerance
%   for convergence of the Forward-Backward sub-iterations (default: 10^-9)
%
% %Output: structure with fields:
% - lambda: the final lambda.
% - A: the mixing matrix.
% - S: the sparse sources matrix.
% With the approximation Y \approx A*S.

function algorithm = ortho_wave_nGMCA_alg(parameters)

if nargin < 1
    parameters = [];
end

algorithm.parameters = parameters;

algorithm.initialize = @m_algo_initialize;
algorithm.iterate = @m_algo_iterate;
algorithm.terminate = @m_algo_terminate;

algorithm.name = 'syn. nGMCA (ortho. wav.)';

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
name = 'ortho_wave_nGMCA_alg';
%main loop
parameters = addMissingParameter(parameters, 'verbose', 0, name);
verbose = parameters.verbose;
parameters = addMissingParameter(parameters, 'phaseRatio',0.85, name, verbose);
parameters = addMissingParameter(parameters, 'MaximumIteration', parameters.rank * 10, name, verbose);

%parameters for A
parameters = addMissingParameter(parameters, 'A.MaximumIteration', 80, name, verbose);
parameters = addMissingParameter(parameters, 'A.RelativeDifferenceTolerance', 0.00001, name, verbose);
%parameters for S
parameters = addMissingParameter(parameters, 'S.tau_MAD', 1, name, verbose);
parameters = addMissingParameter(parameters, 'S.directSparsity', 1, name, verbose);
parameters = addMissingParameter(parameters, 'S.MaximumIteration', 50, name, verbose);
parameters = addMissingParameter(parameters, 'S.qmf', MakeONFilter('Symmlet',4), name, verbose);
parameters = addMissingParameter(parameters, 'S.L', 3, name, verbose);

%% wavelets
if ~exist('FWT_PO', 'file')
    error('ortho_wave_nGMCA_alg requires the WaveLab toolbox to work')
end

n=size(parameters.Y,2);
parameters.S.coarsestLevel=log2(n)-parameters.S.L;
parameters.W = LinearOperator(@(x) FWT_PO(x', parameters.S.coarsestLevel, parameters.S.qmf)',...
    @(x) IWT_PO(x', parameters.S.coarsestLevel, parameters.S.qmf)');


%% random initialization if not provided
if ~isfield(data, 'A')
    data = NMFwarmInit(parameters.Y, parameters.rank);
end




%% first threshold
data = NMFnormalization(data, 'A');
parameters.S.sparsityPattern = ones(1, size(data.S, 2));
if ~parameters.S.directSparsity
    parameters.S.sparsityPattern(: , 1 : (2^parameters.S.coarsestLevel)) = 0;
end
data.lambda = max(al(abs(parameters.W * (data.A' *...
    (data.A * data.S - parameters.Y))))) * ones(size(data.S, 1), 1);
data.S = 0 * data.S;


%% other variables
parameters.refinement_beginning = floor(parameters.phaseRatio * parameters.MaximumIteration);
data.doNotCheckConvergence = 1;

end


%% ----------------------------------------------------------------------%%

function [data, relativeDifference] = m_algo_iterate(data, parameters)
% iterate and update the data

%% update S
AtA = data.A' * data.A;
AtYw = parameters.W * (data.A' * parameters.Y);
Sw = nnWaveSynthesisUpdate(parameters. W * data.S,...
    AtA, AtYw, data.lambda * parameters.S.sparsityPattern,...
    parameters.W, parameters.S);
clear AtA AtYw
data.S = max(0, parameters.W' * Sw);


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

 
 
%% relative difference
relativeDifference = NaN;% save some memory...

%% Update thresholds
data = NMFnormalization(data, 'A');
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
function data = m_thresholdManagement(data, parameters)

 %% compute noise level on the gradient
if isfield(parameters, 'lambdaInf') %if predefine final lambda (not advised)
    sigma_grad = param.lambdaInf / param.S.tau_MAD;
else
    %else, online estimation of the noise on the gradient,
    %based on the MAD estimator of the residual
    
    data.sigma_res = dimMADstd(...
    parameters.W * (data.A * data.S - parameters.Y), 2);
    sigma_grad = sqrt(sum(bsxfun(@times, data.sigma_res.^2, (data.A.^2))))';
    
    % % one may directly compute the residual on the gradient as below
    % % but this is even more biased...
    % sigma_grad = computeWaveScaleFun(...
    % data.A' * (data.A * data.S - parameters.Y),...
    % @(x) parameters.W * x, @(x) dimMADstd(al(x)), parameters.groups);
end


%% decreasing lambda

if data.iteration <= parameters.refinement_beginning
    
    %linear decrease to tau_MAD * sigma_grad when reaching the refinement steps
    data.lambda = max(parameters.S.tau_MAD * sigma_grad, data.lambda ...
        - 1 / (parameters.refinement_beginning + 1 - data.iteration)...
    * (data.lambda - parameters.S.tau_MAD * sigma_grad));
    
    data.doNotCheckConvergence = 1;
else %during refinement, do not modify lambda
    data.doNotCheckConvergence = 0;
end

end


