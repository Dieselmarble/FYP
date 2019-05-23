% nGMCA_oracle_alg.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : J�r�my Rapin (jeremy.rapin@cea.fr)
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
%
% THIS IS AN ORACLE!
% It solves argmin_(S>=0) 1/2*||Y-A_ref S||_2^2+lambda*||S||_1
% and argmin_(A>=0) 1/2*||Y-A S_ref||_2^2
% using iterative soft-thresholding.
% Parameter lambda is set to tau_MAD*sigma_noise, where tau_MAD is a
% constant (preferably in [1,3]) and sigma_noise is estimated using the
% reference data.
%
% 
%
% This function creates a structure triggered by the function ApplyAlgorithm.
% Inputs can be provided through a parameter structure at the creation of the
% algorithm (calling this function) or when using it (calling ApplyAlgorithm).
%
% Required fields of the parameter structure (when calling ApplyAlgorithm if
% not set before):
% - Y: data matrix
% - rank: number of sources
%
% Optional fields:
% - verbose: display information (default:1) or not (0).
% - MaximumFBIteration: maximum number of iteration of the Forward-Backward
%   subroutine (default: 80).
% - FBRelativeDifferenceTolerance: relative difference tolerance
%   for convergence of the Forward-Backward sub-iterations
%   (default: 0.00001).
% - tau_MAD: constant coefficient for the final threshold
%   computation (default: 1).
%
% %Output: structure with fields:
% - lambda: the final lambda.
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
% algo=nGMCA_oracle_alg();
% result = ApplyAlgorithmn(algo,parameters); %perform decomposition
% 
% plot(result.S');
%
%
function algorithm = nGMCA_oracle_alg(parameters)
    if nargin < 1
        parameters = [];
    end
    algorithm.parameters = parameters;

    algorithm.initialize = @m_algo_initialize;
    algorithm.iterate = @m_algo_iterate;
    algorithm.terminate = @m_algo_terminate;

    algorithm.name = 'nGMCA (oracle)';
end



function [data, parameters] = m_algo_initialize(data, parameters)
% initialization of the data structure (and optionnally the function value)

%% check for required parameters
if ~isfield(parameters, 'Y')
    error('In nGMCA_alg, parameter field Y (input data) is required.\n')
end
if ~isfield(parameters, 'rank')
    error('In nGMCA_alg, parameter field rank (number of sources) is required.\n')
end

name = 'nGMCA_oracle_alg';
parameters = addMissingParameter(parameters, 'verbose', 0, name);

parameters = addMissingParameter(parameters, 'S.tau_MAD', 1, name);
parameters = addMissingParameter(parameters, 'S.MaximumIteration', 80, name);
parameters = addMissingParameter(parameters, 'S.RelativeDifferenceTolerance', 0.00001, name);
parameters = addMissingParameter(parameters, 'A.MaximumIteration', 80, name);
parameters = addMissingParameter(parameters, 'A.RelativeDifferenceTolerance', 0.00001, name);


parameters.MaximumIteration = 1;
if isfield(parameters, 'display')
    parameters = rmfield(parameters, 'display');
end
if isfield(parameters, 'recording')
    parameters = rmfield(parameters, 'recording');
end
%% initialization
if ~isfield(parameters, 'reference')
   error('Oracles require reference A and S in parameters.reference.\n'); 
end

parameters.reference = NMFnormalization(parameters.reference, 'A');
data.A = parameters.reference.A;
data.S = parameters.reference.S;

end






function [data, relativeDifference] = m_algo_iterate(data, parameters)
% iterate and update the data (and optionnally the function value)



    sig = 1.4826 * mad(al(parameters.Y -...
        parameters.reference.A * parameters.reference.S)', 1, 2);
    data.lambda = parameters.S.tau_MAD * sig;
    
    reference = NMFnormalization(parameters.reference, 'A');
    % inversion for S
    data.S = nonnegativeSparseUpdate(data.S, reference.A' * reference.A,...
        reference.A' * parameters.Y, data.lambda, parameters.S);
    % inversion for A
    data.A = nonnegativeSparseUpdate(data.A', reference.S * reference.S',...
         reference.S * parameters.Y', 0, parameters.A)';
    
    relativeDifference = 0;
end


function result = m_algo_terminate(data)
% output the result of the algorithm

result.A = data.A;
result.S = data.S;
result.lambda = data.lambda;

end


