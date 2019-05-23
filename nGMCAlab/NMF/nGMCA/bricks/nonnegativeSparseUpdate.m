% nonnegativeSparseUpdate.m - This file is part of nGMCALab.
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
% S = nonnegativeSparseUpdate(S0, AtA, AtY, lambda, parameters)
% Solves problem
% argmin_(S >= 0) 1 / 2 * ||Y - A * S||_2^2 + ||lambda .* S||_1
% using FISTA (Beck & Teboulle 2009)
% Inputs:
% - S0: starting point for S
% - AtA: matrix A' * A
% - AtY: matrix A' * Y
% - lambda: sparsity parameters
% Optional fields for argument "parameter":
%    - MaximumIteration: maximum number of iteration of the algorithm
%    (default: 100).
%    - RelativeDifferenceTolerance: relative difference tolerance (default:
%    0.00001).
%    - normConstrained: 1 to force the norm of the rows of S to be
%    non-increasing, 0 otherwise (default: 0).
%    - normConstraint: specific norm to be enforced for the row norms
%    (default: not used).
%    - reweightedL1: use a reweighting scheme if not 0, using an estimate of noise and
%    the pseudo inverse such that
%    lambda = lambda ./ (1 + (abs(Sinv) ./ noiseEstimate).^2);
%    The noise estimate is computed using the MAD estimator, and multiplied by
%    the value of reweighted L1 (default: 0, typically use 3 otherwise)
% Output:
% S: solution of the problem.
function S = nonnegativeSparseUpdate(S0, AtA, AtY, lambda, parameters)

if nargin < 5
    parameters = [];
end
if isfield(parameters, 'normConstraint')
    parameters.normConstrained = 1;
end
name = 'nonnegativeSparseUpdate';
parameters = addMissingParameter(parameters, 'MaximumIteration', 100, name);
parameters = addMissingParameter(parameters, 'RelativeDifferenceTolerance', 0.00001, name);
parameters = addMissingParameter(parameters, 'normConstrained', 0, name);
parameters = addMissingParameter(parameters, 'reweightedL1', 0, name);

%experimental 
if parameters.reweightedL1 && max(al(lambda)) > 0
   Sinv = AtA \ AtY;
   noiseEst = (parameters.reweightedL1 * dimMADstd(Sinv, 2)) * ones(1, size(Sinv, 2));
   %lambda = lambda .* exp(-abs(S0 ./ noiseEst));
   lambda = lambda ./ (1 + (abs(Sinv) ./ noiseEst).^2);
   %k = 2; plot([S0(k, :); 3 * noiseEst(k, :)]','LineWidth', 2); legend('signal', 'noise', 0);
   clear Sinv noiseEst
end


%precomputation
L = max(eig(AtA));

%initializations
S = S0;
prev_S = S;
t = 1;

%operators
gradient = @(s) AtA * s - AtY;
proximal = @(x, threshold) max(x - threshold, 0);
if parameters.normConstrained
    if isfield(parameters, 'normConstraint')
        norm_limits = parameters.normConstraint;
    else
        norm_limits = dimNorm(S0, 2);
    end
    proximal = @(x, threshold) normProjection(max(x - threshold, 0), norm_limits);
end

for k = 1 : parameters.MaximumIteration
    t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
    w = (t - 1) / t_next;
    t = t_next;
    
    R = (1 + w) * S - w * prev_S;
    prev_S = S;
    S = proximal(R - gradient(R) / L, lambda / L);
    
    if norm(al(prev_S - S), 'fro') / norm(al(S), 'fro') <...
            parameters.RelativeDifferenceTolerance
        break
    end
end

end

