% KLmultNMF_alg.m - This file is part of nGMCALab.
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
%
% This function creates a structure triggered by the function ApplyAlgorithm.
% Inputs can be provided through a parameter structure at the creation of the
% algorithm (calling this function) or when using it (calling
% ApplyAlgorithm).
%
% algorithm = KLmultNMF_alg(parameters)
% Implementation of the multiplicative update algorithm with a
% Kullback-Leibler divergence, proposed in Lee & Seung 2001.

function algorithm = KLmultNMF_alg(parameters)
    if nargin < 1
        parameters = [];
    end
    algorithm.parameters = parameters;

    algorithm.initialize = @m_algo_initialize;
    algorithm.iterate = @m_algo_iterate;
    algorithm.terminate = @m_algo_terminate;

    algorithm.name = 'mult. upd. (KL)';
end



function [data, relativeDifference] = m_algo_iterate(data, parameters)
    %prepare variables
    prevAS = [al(data.A); al(data.S)];
    [m, n] = size(parameters.Y);
    
    %% updates
    data.A = data.A .* ((max(0, parameters.Y) ./ max(eps, data.A * data.S)) * data.S')...
        ./ (ones(m, 1) * max(eps, sum(data.S, 2))');
    data.S = data.S .* (data. A' * (max(0, parameters.Y) ./ max(eps, data.A * data.S)))...
        ./ (max(eps, sum(data.A, 1))' * ones(1, n));
    
    
    %% relative difference
    relativeDifference = norm([al(data.A); al(data.S)] - prevAS) / norm(prevAS);
end

%%
%% NO MODIFICATION REQUIRED BELOW


function [data, parameters] = m_algo_initialize(data, parameters)
% initialization of the data structure (and optionnally the function value)

%% check for required parameters
name = 'KLmultNMF_alg';
if ~isfield(parameters, 'Y')
    error('In %s, parameter field Y (input data) is required.\n', name)
end
if ~isfield(parameters, 'rank')
    error('In %s, parameter field rank (number of sources) is required.\n', name)
end

if ~(isfield(data, 'A') && isfield(data, 'S'))
    r = parameters.rank;
    Y = parameters.Y;
    data = NMFwarmInit(Y, r);
end

% slow algorithm, provide some more iterations
parameters.MaximumIteration = 5 * parameters.MaximumIteration;

end



function result = m_algo_terminate(data)
% output the result of the algorithm
result.A = data.A;
result.S = data.S;

end





