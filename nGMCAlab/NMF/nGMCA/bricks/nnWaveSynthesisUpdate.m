% nnWaveSynthesisUpdate.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
% Created on 7/7/2014, last modified on 16/7/2014
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software. You can use,
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info".
% 
% As a counterpart to the access to the source code and rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty and the software's author,  the holder of the
% economic rights,  and the successive licensors have only limited
% liability.
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean that it is complicated to manipulate,  and that also
% therefore means that it is reserved for developers and experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or
% data to be ensured and,  more generally, to use and operate it in the
% same conditions as regards security.
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.


% Solves argmin_{Sw >= 0} ||Yw - A * Sw * W||_2^2 + ||lambda . * Sw||_1
% using Generalized Forward-Backward algorithm
% Inputs:
% - Sw0: starting point for Sw.
% - AtA: matrix A' * A.
% - AtYw: matrix A' * Yw.
% - lambda: sparsity parameters.
% - W: a tight frame transform of class LinearOperator.
% Optional fields for argument "parameter":
% - MaximumIteration: maximum number of iteration of the algorithm
% (default: 300).
% - RelativeDifferenceTolerance: relative difference tolerance
%   for convergence  (default: 10^-9)
% - display: display function taking "data" as parameter, if one
%  wants to visualize the solutions at each iteration.
%  - record: structure of functions taking "data" as parameter and
%  returning scalar or structure of scalar. They will be
%  recorded in the structure data.recording under the name of the
%  function of its output is a scalar, or of the field if the
%  output is a structure.

function Sw = nnWaveSynthesisUpdate(Sw0, AtA, AtYw, lambda, W, parameters)
% BEWARE: W must be a tight frame !!!!
if nargin < 6
    parameters = [];
end

if ~isfield(parameters, 'MaximumIteration')
    parameters.MaximumIteration = 300;
end
if ~isfield(parameters, 'RelativeDifferenceTolerance')
    parameters.RelativeDifferenceTolerance = 0.000000001;
end

innerParam.gradient = @(x) W * (AtA * (W' * x)) - AtYw;
% can be simplified to @(x) AtA * x - AtYw; in the orthonormal case
innerParam.proximals{1} = @(x, threshold) softThresholding(x, threshold);
innerParam.lambdas{1} = lambda;
innerParam.proximals{2} = @(x, lambda) waveNonnegativity(x, W);
innerParam.lambdas{2} = 0;

innerParam.initialization.x = Sw0;
innerParam.RelativeDifferenceTolerance = parameters.RelativeDifferenceTolerance;
innerParam.MaximumIteration = parameters.MaximumIteration;
innerParam.LipschitzConstant = max(eig(AtA));

if isfield(parameters, 'display')
    innerParam.display = parameters.display;
end
if isfield(parameters, 'recording')
    innerParam.recording = parameters.recording;
end

res = ApplyAlgorithm(GeneralizedForwardBackward_alg, innerParam);
Sw = res.x;

end

function yw = waveNonnegativity(xw, W)
% this is where the tight frame requirement is used (straightforward proximal
% operator)
    z = max(0, -(W' * xw));
    yw = xw + (W * z);
% can be simplified to yw = W * max(0, W' * xw); in the orthonormal case
end


