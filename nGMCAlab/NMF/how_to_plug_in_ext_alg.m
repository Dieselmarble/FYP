% how_to_plug_in_ext_alg.m - This file is part of nGMCALab.
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
% algorithm = how_to_plug_in_alg(parameters)
% Function showing how to create a structure for use with ApplyAlgorithm
% from an external algorithm. Since it bypasses all the features of
% ApplyAlgorithm, the main interest is to be able to use any external
% algorithm with the performBenchmark tool.


function algorithm = how_to_plug_in_ext_alg(parameters)
    if nargin < 1
        parameters = [];
    end
    algorithm.parameters = parameters;

    algorithm.initialize = @m_algo_initialize;
    algorithm.iterate = @m_algo_iterate;
    algorithm.terminate = @m_algo_terminate;

    algorithm.name = 'modify_name';%TO MODIFY
end



function [data, relativeDifference] = m_algo_iterate(data, parameters)
    relativeDifference = 0;
    
    r = parameters.rank;
    Y = parameters.Y;
        
    %if parameters.reference is specified (performBenchmark.m specifies
    %it), it can be used here for parameter selection in case of Oracles. 
    
    [A, S] = nGMCA(Y, r);% MODIFY TO ANY OTHER ALGORITHM
    
    data.A = A;
	data.S = S;


end

%%
%% NO MODIFICATION REQUIRED BELOW


function [data, parameters] = m_algo_initialize(data, parameters)
% initialization of the data structure (and optionnally the function value)

%% check for required parameters
name = 'how_to_plugin_in_ext_alg';
if ~isfield(parameters, 'Y')
    error('In %s, parameter field Y (input data) is required.\n', name)
end
if ~isfield(parameters, 'rank')
    error('In %s, parameter field rank (number of sources) is required.\n', name)
end

%% BYPASS APPLYALGORITHM FEATURES
parameters.MaximumIteration = 1; 
if isfield(parameters, 'display')
    parameters = rmfield(parameters, 'display');
end
if isfield(parameters, 'recording')
    parameters = rmfield(parameters, 'recording');
end
%make sure there is only one iteration
parameters.MaximumIteration = 1;



end



function result = m_algo_terminate(data)
% output the result of the algorithm
result.A = data.A;
result.S = data.S;

end





