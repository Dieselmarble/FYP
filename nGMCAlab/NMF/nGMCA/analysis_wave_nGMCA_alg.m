% analysis_wave_nGMCA_alg.m - This file is part of nGMCALab.
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
% non-negative Generalized Morphological Component Analysis with analysis
% sparse regularization in the wavelet domain
%
% Aims at solving :
% argmin_(A >= 0, S >= 0) 1 / 2 * || Y - A * S||_2^2 + lambda * ||S * W^T||_1
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
% - S.sourcesShape: when dealing with 2D images, shape of the images (ex: [128, 128])
% - S.tau_MAD: constant coefficient for the final threshold computation (default: 2).
% - S.directSparsity: penalize the coarse scale or not (default: 1)
% - S.reweightedL1: use a reweighting scheme if non-null (default: 0)
% - S.qmf: wavelet coefficients (use MakeONFilter from the Wavelab toolbox,
% default: Symmlet, 4)
% - S.L: number of scales of the wavelets (default: 3)
% - S.isometric: use isometric wavelet transform (default: 0)
% - S.uniformFirstThreshold: use the same threshold for each source at the
% first iteration (default: 1). If the sources have very different
% dynamics, better set to 0.
% Optional subfields for A:
% - A.MaximumIteration: number of iterations of the update of S  (default: 80)
% - A.RelativeDifferenceTolerance: relative difference tolerance
%   for convergence of the Forward-Backward sub-iterations (default: 10^-9)
% - A.noL2CoarseScale: do not use the wavelet coarse scale of S and Y for
% the update of A (default: 0)
%
% %Output: structure with fields:
% - lambda: the final lambda.
% - A: the mixing matrix.
% - S: the sparse sources matrix.
% With the approximation Y \approx A*S.

function algorithm = analysis_wave_nGMCA_alg(parameters)

if nargin < 1
    parameters = [];
end

algorithm.parameters = parameters;

algorithm.initialize = @m_algo_initialize;
algorithm.iterate = @m_algo_iterate;
algorithm.terminate = @m_algo_terminate;

algorithm.name = 'ana. nGMCA (red. wav.)';

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
name = 'analysis_wave_nGMCA_alg';
%main loop
parameters = addMissingParameter(parameters, 'verbose', 0, name);
verbose = parameters.verbose;
parameters = addMissingParameter(parameters, 'phaseRatio',0.85, name, verbose);
parameters = addMissingParameter(parameters, 'MaximumIteration', 500, name, verbose);

%parameters for A
parameters = addMissingParameter(parameters, 'A.MaximumIteration', 80, name, verbose);
parameters = addMissingParameter(parameters, 'A.RelativeDifferenceTolerance', 0.000000001, name, verbose);
parameters = addMissingParameter(parameters, 'A.noL2CoarseScale', 0, name, verbose);

%parameters for S
parameters = addMissingParameter(parameters, 'S.tau_MAD', 2, name, verbose);
parameters = addMissingParameter(parameters, 'S.sourcesShape', size(parameters.Y, 2), name, 0);
parameters = addMissingParameter(parameters, 'S.directSparsity', 1, name, verbose);
parameters = addMissingParameter(parameters, 'S.reweightedL1', 0, name, verbose);
parameters = addMissingParameter(parameters, 'S.MaximumIteration', 24, name, verbose);
%filters must be wavelet filters. Use the MakeONFilter of the WaveLab toolbox
parameters = addMissingParameter(parameters, 'S.qmf', MakeONFilter('Symmlet',4), name, verbose);
parameters = addMissingParameter(parameters, 'S.L', 3, name, verbose);
parameters = addMissingParameter(parameters, 'S.isometric', 0, name, verbose);
parameters = addMissingParameter(parameters, 'S.uniformFirstThreshold', 1, name, verbose);
parameters.S.nonNegative = 1; %for the redWave toolbox
parameters.S.groups = []; %to create specific wavelet groups (depleted)

%% random initialization if not provided
if ~isfield(data, 'A')
    data = NMFwarmInit(parameters.Y, parameters.rank);
end

%% prepare wavelet operator
if isscalar(parameters.S.sourcesShape)
    parameters.S.dimensions = 2;
    parameters.shape2line = @(S) S;
    parameters.line2shape = @(S) S;
    parameters.W = LinearOperator(@(x) redWave(x, 1, parameters.S),...
        @(x) redWave(x, -1, parameters.S));
else
    parameters.S.dimensions = 2 : (length(parameters.S.sourcesShape) + 1);
    
    parameters.shape2line = @(S) reshape(S,...
        [size(S, 1), size(S, 2) * size(S, 3)]);
    parameters.line2shape = @(S) reshape(S,...
        [size(S, 1), parameters.S.sourcesShape]);
    
    parameters.W = LinearOperator(@(x) redWave(parameters.line2shape(x),...
        1, parameters.S),...
        @(x) parameters.shape2line(redWave(x, -1, parameters.S)));
    
end


%% first threshold
data = NMFnormalization(data, 'A');

% compute maximum on each scale
%%data.lambda = computeWaveScaleFun(...
%%    data.A' * (data.A * data.S - parameters.Y),...
%%    @(x) parameters.W * x, @(x) max(al(abs(x))), parameters.S.groups);
data.lambda = computeWaveScaleFun([size(data.S, 1), parameters.S.sourcesShape],...
    parameters.W * (data.A' * (data.A * data.S - parameters.Y)),...
    @(x) max(al(abs(x))), parameters.S.groups);

% uniformize among sources
% this may not be helpful if the sources have very different dynamics
if parameters.S.uniformFirstThreshold
data.lambda = bsxfun(@times,...
    ones(parameters.rank, 1), max(data.lambda, [], 1));
end
% remove the threshold on the coarse scale if not sparse in the direct
% domain
if isempty(parameters.S.groups)
    data.lambda(:, 1, 1) = (parameters.S.directSparsity ~= 0) * data.lambda(:, 1, 1);
end
data.S = 0 * data.S;


%% other variables
parameters.refinement_beginning = floor(parameters.phaseRatio * parameters.MaximumIteration);
data.doNotCheckConvergence = 1;


%% reweighting
if parameters.A.noL2CoarseScale
    param = parameters.S;
    param.isometric = 1;
    parameters.A.Wn = LinearOperator(@(x) redWave(parameters.line2shape(x),...
        1, param),...
        @(x) parameters.shape2line(redWave(x, -1, param)));
    parameters.A.Ywn = parameters.shape2line(parameters.A.Wn * parameters.Y);
    
    
    if isscalar(parameters.S.sourcesShape)
        parameters.A.mask = ones(1, parameters.S.sourcesShape * (parameters.S.L + 1));
        parameters.A.mask(1, 1 : size(data.S, 2)) = 0;
    else
        sizes = [1, parameters.S.sourcesShape];
        sizes(parameters.S.dimensions) = 2 * sizes(parameters.S.dimensions);
        sizes(parameters.S.dimensions(1)) = parameters.S.L * sizes(parameters.S.dimensions(1));
        parameters.A.mask = ones(sizes);
        n1 = parameters.S.sourcesShape(1);
        n2 = parameters.S.sourcesShape(2);
        d1 = parameters.S.dimensions(1);
        for k = 1 : parameters.S.L
            parameters.A.mask(1, (1 : n1) + 2 * (d1 == 2) * (k - 1) * n1,...
                (1 : n2) + 2 * (d1 == 3) * (k - 1) * n2) = 0;
        end
        parameters.A.mask = parameters.shape2line(parameters.A.mask);
    end
    parameters.A.Ywn = bsxfun(@times, parameters.A.Ywn, parameters.A.mask);
        
end

end


%% ----------------------------------------------------------------------%%

function [data, relativeDifference] = m_algo_iterate(data, parameters)
% iterate and update the data


%% threshold
Lambdas = generateFullGroupMatrix(data.lambda, [1, parameters.S.sourcesShape], parameters.S.groups);
if parameters.S.reweightedL1
    S_inv = data.A \ parameters.Y;
    data.S_inv = S_inv;
    S_w = parameters.W * S_inv;
    
    sig_sources = computeWaveScaleFun([size(data.S, 1), parameters.S.sourcesShape],...
        parameters.W * S_inv, @(x) dimMADstd(al(x)), parameters.S.groups);
    Sigmas = generateFullGroupMatrix( sig_sources, [1, parameters.S.sourcesShape], parameters.S.groups);
    
    Lambdas = Lambdas ./ (1 + (abs(S_w) ./ Sigmas).^2);
end
data.Lambdas = Lambdas;

%% update S
AtA = data.A' * data.A;
AtY = parameters.line2shape(data.A' * parameters.Y);
S = parameters.line2shape(data.S);
S = redWaveInversion(S, AtA, AtY,...
    Lambdas, parameters.S);
data.S = parameters.shape2line(S);
clear AtA AtY S


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

if ~parameters.A.noL2CoarseScale
    SSt = data.S * data.S';
    SYt = data.S * parameters.Y';
else
    Swn = parameters.shape2line(parameters.A.Wn * data.S);
    Swn = bsxfun(@times, Swn, parameters.A.mask);
    SSt = Swn * Swn';
    SYt = Swn * parameters.A.Ywn';
end

data.A = nonnegativeSparseUpdate(data.A', SSt,...
    SYt, 0, parameters.A)';

%% relative difference
data = NMFnormalization(data, 'A');
relativeDifference = NaN;% save some memory...


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
function data = m_thresholdManagement(data, parameters)
sigma_update_ratio = 0.05;

if isfield(parameters, 'addStd') %predefined level of noise
    if isscalar(parameters.addStd)
        data.sigma_res = parameters.addStd * ...
            ones(size(data.A, 1), size(data.lambda, 2), size(data.lambda, 3));
    else
        data.sigma_res = parameters.addStd;
    end
else
%% compute noise level on the gradient

% online estimation of the noise on the gradient,
% based on the MAD estimator of the residual
if data.iteration == 1 % not necessary to update all the residues at once
    ind = 1 : size(data.A, 1);
    data.sigma_res = zeros(size(data.A, 1), size(data.lambda, 2),...
        size(data.lambda, 3));
    for k = ind
        %data.sigma_res(k, :, :) = computeWaveScaleFun(...
        %data.A(k, :) * data.S - parameters.Y(k, :),...
        %@(x) parameters.W * x, @(x) dimMADstd(al(x)), parameters.S.groups);
        data.sigma_res(k, :, :) = computeWaveScaleFun([1, parameters.S.sourcesShape],...
        parameters.W * (data.A(k, :) * data.S - parameters.Y(k, :)),...
        @(x) dimMADstd(al(x)), parameters.S.groups);
    end
else
    ind = randperm(size(data.A, 1));
    ind = ind(1 : floor(size(data.A, 1) * sigma_update_ratio + 0.99));
    %data.sigma_res(ind, :, :) = computeWaveScaleFun(...
    %    data.A(ind, :) * data.S - parameters.Y(ind, :),...
    %    @(x) parameters.W * x, @(x) dimMADstd(al(x)), parameters.S.groups);
    data.sigma_res(ind, :, :) = computeWaveScaleFun([numel(ind), parameters.S.sourcesShape],...
    parameters.W * (data.A(ind, :) * data.S - parameters.Y(ind, :)),...
    @(x) dimMADstd(al(x)), parameters.S.groups);

end
%     data.sigma_res(ind, :, :) = computeWaveScaleFun(...
%     data.A(ind, :) * data.S - parameters.Y(ind, :),...
%     @(x) parameters.W * x, @(x) dimMADstd(al(x)), parameters.S.groups);
end

sigma_grad = data.sigma_res(1 : size(data.S, 1), :, :);
for k = 1 : size(data.S, 1)
    sigma_grad(k, :, :) = sqrt(sum(bsxfun(@times, data.A(:, k).^2, data.sigma_res.^2)));
end

% % one may directly compute the residual on the gradient as below
% % but this is even more biased...
% sigma_grad = computeWaveScaleFun(...
% data.A' * (data.A * data.S - parameters.Y),...
% @(x) parameters.W * x, @(x) dimMADstd(al(x)), parameters.groups);

%% decreasing lambda

if data.iteration <= parameters.refinement_beginning
    
    %linear decrease to tau_MAD * sigma_grad when reaching the refinement steps
    data.lambda = max(parameters.S.tau_MAD * sigma_grad, data.lambda ...
        - 1 / (parameters.refinement_beginning + 1 - data.iteration)...
        * (data.lambda - parameters.S.tau_MAD * sigma_grad));
    data.lambda(:, 1, 1) = (parameters.S.directSparsity ~= 0) * data.lambda(:, 1, 1);
    
    data.doNotCheckConvergence = 1;
else %during refinement, do not modify lambda
    data.doNotCheckConvergence = 0;
end

end



%%
%%



