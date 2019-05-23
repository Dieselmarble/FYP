% GeneralizedForwardBackward_alg.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
% Created on 13/6/2014, last modified on 16/7/2014
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

% function algorithm = GeneralizedForwardBackward_alg(parameters)
% This function creates a structure triggered by the function ApplyAlgorithm.
% Inputs can be provided through a parameter structure at the creation of the
% algorithm (calling this function) or when using it (calling ApplyAlgorithm).
%
% The Generalized Forward-Backward algorithm (see Raguet et al. 2001) solves
% min_x f(x)+\sum_i \lambda_i g_i(x) with f convex differentiable with Lipschitz gradient and
% g_i proper convex and lower semi-continuous functions.
% It requires the following fields in the parameter structure (when calling
% ApplyAlgorithm if not set before):
% - proximals: cell list of the proximal operators of g_i under the form @(x,lambda) prox.
% - lambdas: cell list of the lambda parameters in the above problem
% - gradient: gradient of f under the form @(x) grad.
% - LipschitzConstant: Lipschitz constant of the gradient of f.
% - initialization.x: initialization of the algorithm.


function algorithm  =  GeneralizedForwardBackward_alg(parameters)

% if no parameter is given,  used default ones
parameters.name = 'GeneralizedForwardBackward';
algorithm.parameters = parameters;

%%
algorithm.initialize = @m_algo_initialize;
algorithm.iterate = @m_algo_iterate;
algorithm.terminate = @m_algo_terminate;
algorithm.name = 'Generalized Forward-Backward';

end



%% ----------------------------------------------------------------------%%

function [data, parameters] = m_algo_initialize(data, parameters)
% initialization of the algorithm


%% check for required parameters
if ~isfield(parameters, 'proximals')
    error('In GeneralizedForwardBackward_alg: You need to input the proximal operators cell in the parameter structure.\n')
end
if ~isfield(parameters, 'gradient')
    error('In GeneralizedForwardBackward_alg: You need to input the gradient function in the parameter structure.\n')
end
if ~isfield(parameters, 'LipschitzConstant')
    error('In GeneralizedForwardBackward_alg: You need to input LipschitzConstant for the gradient in the parameter structure (2*lambda_max).\n')
end
if ~isfield(parameters, 'lambdas')
    error('In GeneralizedForwardBackward_alg: You need to input lambdas cell in the parameter structure (2*lambda_max).\n')
end
if ~isfield(data, 'x')
    error('In GeneralizedForwardBackward_alg: You need to input initialization.x in the parameter structure.\n')
end

%% initialize automatically missing general parameters
name = 'GeneralizedForwardBackward_alg';
parameters = addMissingParameter(parameters, 'verbose', 0, name, 0);
parameters = addMissingParameter(parameters, 'MaximumIteration', 1000, name, parameters.verbose);
parameters = addMissingParameter(parameters, 'MaximumIteration', 10000, name, parameters.verbose);
parameters = addMissingParameter(parameters, 'RelativeDifferenceTolerance', 10^-6, name, parameters.verbose);


%% initialization
parameters.NumProximals = length(parameters.proximals);
data.t = 1;
for k = 1 : parameters.NumProximals
   data.z{k} = data.x; 
end

end


%% ----------------------------------------------------------------------%%

function [data, relativeDifference] = m_algo_iterate(data, parameters)
% iterate and update the data (and optionnally the function value)


%%
omega_i = 1 / parameters.NumProximals;
beta = 1 / parameters.LipschitzConstant;
gamma_t = beta;
mu_t = 0.9 * min(3 / 2, (1 + 2 * beta / gamma_t) / 2); %lambda_t in the paper
grad = parameters.gradient;


%%
x  =  data.x;
data.x = 0 * data.x;

temp = 2 * x-gamma_t * grad(x);
for i = 1:parameters.NumProximals;

    prox = parameters.proximals{i};
    data.z{i} = data.z{i} + mu_t * (prox(temp - data.z{i}, ...
        parameters.lambdas{i} * gamma_t / omega_i) - x);
    
    %update x
    data.x = data.x + omega_i * data.z{i};
end



%%
relativeDifference = norm(al(data.x - x)) / norm(al(x));


end



%% ----------------------------------------------------------------------%%

function result = m_algo_terminate(data)
% output the result of the algorithm

result.x = data.x;
end



