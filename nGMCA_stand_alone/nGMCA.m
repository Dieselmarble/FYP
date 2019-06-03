% nGMCA.m - This file is part of nGMCALab (stand alone file)
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13
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
% Aims at solving argmin_(A >= 0, S >= 0) 1 / 2 * ||Y - AS||_2^2 + lambda * ||S||_1
% using iterative soft-thresholding.
% Parameter lambda is decreasing during the iterations and set to
% tau_MAD*sigma_noise at the end of the algorithm, where tau_MAD is a
% constant (preferably in [1,3]) and sigma_noise is an online estimate of
% the noise standard deviation.
%
% For more information, this algorithm is described as nGMCA^S in:
% J. Rapin, J.Bobin, A. Larue and J.L. Starck,
% Sparse and Non-negative BSS for Noisy Data,
% IEEE Transactions on Signal Processing, 2013.
% Please use the above reference if using this code in a publication.
% 
%
% %Inputs:
% - Y: data matrix
% - r: number of sources
% - parameters: optional structure for advanced parametring:
%           - verbose: display informations (default:1) or not (0).
%           - MaximumIteration: maximum number of iterations of the main loop
%           (default: 500, should be enough for convergence).
%           - MaximumFBIteration: maximum number of iteration of the Forward-Backward
%             subroutine (default: 80).
%           - FBRelativeDifferenceTolerance: relative difference tolerance
%           for convergence of the Forward-Backward sub-iterations
%           (default: 0.00001).
%           - phaseRatio: transition between decreasing thresholding phase
%           and refinement phase in pourcent of the iterations (default: 0.80)
%           - tau_MAD: constant coefficient for the final threshold
%           computation (default: 1).
%           - lambdaInf: in order to predefine the final lambda (not advised,
%           except for underdetermined settings).
%
% %Outputs:
% - A: the mixing matrix.
% - S: the sparse sources matrix.
% with the approximation Y \approx A*S.
%
% %Example of use:
% m = 200; %number of observations
% n = 200; %number of source samples
% r = 5; %number of sources
% ref_A = rand(m, r); %non-negative mixtures
% ref_S = rand(r, m) .* (rand(r, m) > 0.9); %non-negative and sparse sources
% Y = ref_A * ref_S + 0.1 * randn(m, n); %simulated data
% [A,S] = nGMCA(Y, r); %perform decomposition
% 
% %same but without text
% parameters.verbose = 0;
% [A,S] = nGMCA(Y, r, parameters);
% %more options are available, check help for more information
% 
% %plot
% subplot(2, 1, 1); plot(ref_S' / max(max(ref_S)), 'LineWidth', 2);
% set(gca, 'FontName', 'TimesNewsRoman');
% title('Reference sources', 'FontName', 'TimesNewsRoman');
% subplot(2,1,2); plot(S' / max(max(S)), 'LineWidth', 2);
% title('Estimated sources (up to scale and permutation invariances)', 'FontName', 'TimesNewsRoman');
% set(gca, 'FontName', 'TimesNewsRoman');
% % colors do not match because the sources are not reordered.
%
%
% N.B.: - This is the main algorithm of the paper, which is prefered to the
%       hard-thresholding and naive approaches.
%
function [A,S] = nGMCA(Y, r, parameters)


%% Check and initialize parameters
parameters.name = 'nGMCA';
parameters = m_automaticSettings('verbose', 1, parameters);
parameters = m_automaticSettings('MaximumIteration', 500, parameters);
parameters = m_automaticSettings('MaximumFBIteration', 80, parameters);
parameters = m_automaticSettings('FBRelativeDifferenceTolerance', 0.00001, parameters);
parameters = m_automaticSettings('phaseRatio', 0.80, parameters);
parameters = m_automaticSettings('tau_MAD', 1, parameters);
refinement_beginning = floor(parameters.phaseRatio * parameters.MaximumIteration);


%% initialize data
%initialization and warm up iterations
data = m_warmInitNMF(Y, r);

%first threshold
data = m_normalizeA(data, 1);
data.lambda = max(max(abs(data.A' * (data.A * data.S - Y))));


%% iterate

    if parameters.verbose
        str = [];
        fprintf(1, 'Iteration ');
    end    
for iter = 1 : parameters.MaximumIteration
    data.iteration = iter;

    %% update S
    data.S = m_nonnegativeSparseUpdate_S(data.A, data.S, Y, data.lambda, parameters);

    %reinitialize lines of S if need be
    data = m_reinitializeS(data, Y, parameters.verbose);
        
    
    %% update A
    %preconditioning
    data = m_normalizeA(data, 0);
    %A column norms must be non-increasing in order to converge to a
    %stationary point
    %this is only helpful in the refinement steps (since the cost function
    %is then fixed)
    if iter >= refinement_beginning
        parameters.normConstrained = 1;
    end
    %update
    data.A = m_nonnegativeSparseUpdate_A(data.A, data.S, Y, 0, parameters);

    
    %% update lambda
    data = m_normalizeA(data, 1);
    data.lambda = m_updateLambda(data, Y, parameters);
    
    %%
    if parameters.verbose
        fprintf(repmat('\b', 1, length(str)));
        str = [num2str(data.iteration) '/' num2str(parameters.MaximumIteration) '.'];
    	fprintf(1, str);
    end    
end

if parameters.verbose
    fprintf(1, '\n');
end    
%% Output
A = data.A;
S = data.S;

end



%% initialization
function data = m_warmInitNMF(Y, r)
% initialization of the data structure
% alternates between a couple of unconstrained (fast) update and constrained update

    m = size(Y, 1);
    data.A = rand(m, r);
    options.MaximumFBIteration = 50;
    options.FBRelativeDifferenceTolerance = 0;
    
    for k = 1 : 2
        
        %% unconstrained updates
        i = sum(data.A) > 0;
        data.S(i, :) = max(data.A(:, i) \ Y, 0);
        i = sum(data.S, 2) > 0;
        data.A(:, i)  =  max(Y / data.S(i, :), 0);


        %% constrained updates
        data.S = m_nonnegativeSparseUpdate_S(data.A, data.S, Y, 0, options);
        % reinitialize S if need be
        data = m_reinitializeS(data, Y, 0);
        data.A = m_nonnegativeSparseUpdate_A(data.A, data.S, Y, 0, options);
        
    end
end





%% optimization procedures
function S = m_nonnegativeSparseUpdate_S(A, S0, Y, lambda, parameters)
%solves problem argmin_(S >= 0) 1 / 2 * ||Y - A * S||_2^2 + lambda * ||S||_1
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
proximal = @(x,threshold) max(x - threshold, 0);

for k = 1 : parameters.MaximumFBIteration
    t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
    w = (t - 1) / t_next;
    t = t_next;
    
    R = (1 + w) * S - w * prev_S;
    prev_S = S;
    S = proximal(R - gradient(R) / L, lambda / L);

    if norm(al(prev_S - S), 'fro') / norm(al(S), 'fro') < parameters.FBRelativeDifferenceTolerance
        break
    end
end

end

function A = m_nonnegativeSparseUpdate_A(A0, S, Y, lambda, parameters)
%solves problem argmin_(A >= 0) 1 / 2 * ||Y - A * S||_2^2
%with starting point A0
%can optionally constrain the norm of the columns of A to be
%non-increasing

%precomputation
H = S * S';
YSt = Y * S';
L = max(eig(H));

%initializations
A = A0;
prev_A = A;
t = 1;

%operators
gradient = @(a) a * H - YSt;
proximal = @(x,threshold) max(x, 0);
if isfield(parameters, 'normConstrained')
    if parameters.normConstrained
        norm_limits = m_dimNorm(A0, 1);
        proximal = @(x,threshold) m_normProj(max(x, 0), norm_limits);
    end
end


for k = 1 : parameters.MaximumFBIteration
    t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
    w = (t - 1) / t_next;
    t = t_next;
    
    B = (1 + w) * A - w * prev_A;
    prev_A = A;
    A = proximal(B - gradient(B) / L, lambda / L);

    if norm(al(prev_A - A), 'fro') / norm(al(A), 'fro') < parameters.FBRelativeDifferenceTolerance
        break
    end
end

end







%% reinitialization when a source was set to 0

function data = m_reinitializeS(data, Y, verbose)
%straightforward reinitialization of a line of S and a column of A
%by picking one column in the residue

if sum(sum(data.S, 2) == 0) > 0 || sum(sum(isnan(data.S))) > 0
    if verbose
        fprintf(1, 'Reinitialization of %i null line(s) in S.\n', sum(sum(data.S, 2) == 0));
    end
    indices = ((sum(data.S, 2) == 0) + sum(isnan(data.S), 2) > 0) > 0;
    [data.A(:, indices), data.S(indices, :)] = m_fastExtractNMF(Y - data.A * data.S, sum(indices));
end
end


function [A,S] = m_fastExtractNMF(residual, r)

if r > 0
    [m, n] = size(residual);
    A = zeros(m, r);
    S = zeros(r, n);
    for i = 1 : r
        residual = residual .* (residual > 0);
        if sum(residual(:)) ~= 0
            % compute square norm of residual to select maximum one
            res2 = sum(residual .* residual);
            j = find(res2 == max(res2));
            if ~isempty(j)
                j = j(1);
                if res2(j) > 0
                    %normalize maximum residual
                    A(:, i) = residual(:, j) / sqrt(res2(j));
                    %compute scalar product with the rest od the residual and keep only
                    %positive coefficients
                    S(i, :) = A(:, i)' * residual;
                    S(i, :) = max(S(i, :), 0);
                    %compute new residual
                    residual = residual - A(:, i) * S(i, :);
                end
            end
        end
    end
else
    A = [];
    S = [];
end
end


%% Tools
function parameters = m_automaticSettings(name, value, parameters)
%handling of automatic parameters and display

if ~isfield(parameters,name)
    parameters.(name) = value;
    if isfield(parameters, 'verbose')
        if parameters.verbose
            fprintf(1, 'In %s, parameter "%s" automatically set to %f.\n', parameters.name,name,parameters.(name));
        end
    end
end
end

function data = m_normalizeA(data,true)
% normalization of the NMF data
% Inputs:
% - data: structure with fields A and S (size(data.A, 2) == size(data.S, 1))
% - true: if true == 1, A columns are normalized and the scale is on S
%         and conversely for true == 0
% Ourput:
% - data: structure with fields A and S (size(data.A, 2) == size(data.S, 1))

    Ds = m_dimNorm(data.S', 1);
    Da = m_dimNorm(data.A, 1);
    D = Da .* Ds;
    act = D > 0;
    if true == 0      
        data.S(act, :) = bsxfun(@times, data.S(act, :), 1 ./ Ds(act)');
        data.A(:, act) = bsxfun(@times, data.A(:, act), D(act) ./ Da(act));
    else
        data.S(act, :) = bsxfun(@times, data.S(act, :), (D(act) ./ Ds(act))');
        data.A(:, act) = bsxfun(@times, data.A(:, act), 1 ./ Da(act));
    end
    data.S(~act, :) = 0 * data.S(~act, :);
    data.A(:, ~act) = 0 * data.A(:, ~act);
end

function norms = m_dimNorm(x, dim)
%Inputs:
%- x: tensor data
%- dim: a dimension (integer)
%Output:
%- norms: L2 norm of the columns (dim = 1) or lines (dim = 2) of x

norms = sqrt(sum(x .* x, dim));

end



function X = m_normProj(X, norm_limits)
%Inputs:
%- X: matrix
%- norm_limits: maximum norm allowed for the lines (if norm_limits is a 
%  column vector) or columns (if norm_limits is a line vector) of X
%Output:
%- X: matrix after projection on the norm constraints

m = size(X, 1);
if size(norm_limits, 1) == m
    norms = m_dimNorm(X, 2);
    i = norms > 0;
    X(i, :) = bsxfun(@times, X(i, :), min(norms(i), norm_limits(i)) ./ norms(i));
    
else
    norms = m_dimNorm(X, 1);
    i = norms > 0;
    X(:, i) = bsxfun(@times, X(:, i), min(norms(i), norm_limits(i)) ./ norms(i));
end


end




%% Update sparsity parameter
function lambda = m_updateLambda(data, Y, param)
refinement_beginning = floor(param.phaseRatio * param.MaximumIteration);

if refinement_beginning > data.iteration
    
    if isfield(param, 'lambdaInf') %if predefine final lambda (not advised)
        sigma_residue = param.lambdaInf / param.tau_MAD;
    else %else, online estimation based on the MAD estimator of std of the residue
        %(since A columns are normalized to L2 unity)
        sigma_residue = m_dimMADstd(al(Y - data.A * data.S), 1);
    end

    %linear decrease to tau_MAD*sigma_residue when reaching the refinement steps
    lambda = max(param.tau_MAD * sigma_residue, data.lambda ...
        - 1/(refinement_beginning - data.iteration)...
        * (data.lambda - param.tau_MAD * sigma_residue));
else %during refinement, do not modify lambda
    lambda = data.lambda;
end

end

function sigmas = m_dimMADstd(X, dim)
% compute the standard deviation of X based on the MAD estimator
% along dimension dim

medval = median(X, dim);
sigmas = 1.4826 * median(abs(bsxfun(@minus, X, medval)), dim);

end

function y = al(x)
    y = x(:);
end
