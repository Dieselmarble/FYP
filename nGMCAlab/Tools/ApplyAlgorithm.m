% ApplyAlgorithm.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
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
% result = ApplyAlgorithm(optimizationAlgorithm, parameters)
% Inputs:
% - optimizationAlgorithm: algorithm structure with fields initialize, 
% iterate and terminate (see ForwardBackward_struct for an example)
% - parameter: structure containing parameters of the algorithm.
%       Required parameters depend on the type of algorithm (see algorithms
%       code and description to find which are necessary).
%       Common parameters include:
%           - MaximumIteration: maximum number of steps to perform
%           (default: 10^6)
%           - RelativeDifferenceTolerance: tolerance of difference between
%           an iterate and the next (default: 100*eps)
%           - MaximumTime: maximum time for running the algorithm (default:
%           not used)
%           - display: display function taking "data" as parameter, if one
%           wants to visualize the solutions at each iteration.
%           - record: structure of functions taking "data" as parameter and
%           returning scalar or structure of scalar. They will be
%           recorded in the structure data.recording under the name of the
%           function of its output is a scalar, or of the field if the
%           output is a structure.
% Output: 
% - result: structure with fields returned by optimizationAlgorithm.terminate
%       It also contains the structure recording which fields are the values
%       of the recording function at each iterate.
%
%
%

function result = ApplyAlgorithm(optimizationAlgorithm, parameters)


% initialize parameter structure if need be
if nargin < 2
    parameters = [];
end

if ~isfield(optimizationAlgorithm, 'initialize')
    error('Algorithm structure requires a field "initialize".\n')
end
if ~isfield(optimizationAlgorithm, 'iterate')
    error('Algorithm structure requires a field "iterate".\n')
end
if ~isfield(optimizationAlgorithm, 'terminate')
    error('Algorithm structure requires a field "terminate".\n')
end



%for convenient use of the name within the algorithm
if ~isfield(optimizationAlgorithm, 'name')
    error('Algorithm structure requires a field "name".\n')
else
    parameters.name = optimizationAlgorithm.name;
end




%% INITIALIZATION OF THE ALGORITHM
% combine parameters given to the framework and to the algorithm
% priority is given to the framework in case of double affectation
parameters = m_combineFields(parameters, optimizationAlgorithm.parameters);


% data initialization
% use initialization if it exists
if isfield(parameters, 'initialization')
    data = parameters.initialization;
else
    data = [];
end

%algorithm initialization
[data, parameters] = optimizationAlgorithm.initialize(data, parameters);

%complete initialization if need be
parameters = addMissingParameter(parameters, 'RelativeDifferenceTolerance', 100 * eps, 'ApplyAlgorithm');
parameters = addMissingParameter(parameters, 'MaximumIteration', 10^6, ' ApplyAlgorithm');
parameters.MaximumIteration = floor(parameters.MaximumIteration);

%% ITERATIONS

data.recording.relativeDifference = NaN * ones(parameters.MaximumIteration, 1);
keepLooping = true;
data.iteration = 1;
if parameters.MaximumIteration < 1
    keepLooping = false;
end
time_tic = tic;
while keepLooping
    
    
    %% iterate
    [data, relativeDifference] = optimizationAlgorithm.iterate(data, parameters);
    time = toc(time_tic);
    
    % recording relative difference
    data.recording.relativeDifference(data.iteration) = relativeDifference;
    
    
    %% recording values
    if isfield(parameters, 'recording')
        data = m_record(data, parameters);
    end
    
    %% display
    if isfield(parameters, 'display')
        parameters.display(data);
        pause(0.01);
        refresh;
        if ~isnan(relativeDifference)
            fprintf(1, 'Relative difference at iteration %i/%i is %f%%.\n',...
                data.iteration, parameters.MaximumIteration, 100 * relativeDifference);
        else
            fprintf(1, 'Iteration %i/%i\n',...
                data.iteration, parameters.MaximumIteration);
        end
        
        if isfield(data, 'recording')
            if (isfield(data.recording, 'SDR_S') && isfield(data.recording, 'SDR_A'))
                fprintf('SDR_S = %f, SDR_A = %f.\n',...
                    data.recording.SDR_S(data.iteration),...
                    data.recording.SDR_A(data.iteration))
            end
        end
        
    end
    
    keepLooping = m_checkConvergence(data, parameters, relativeDifference, time);
    data.iteration=data.iteration + 1;
end

%% terminate algorithm

if nargin(optimizationAlgorithm.terminate) == 1
    result = optimizationAlgorithm.terminate(data);
else
    result = optimizationAlgorithm.terminate(data, parameters);
end


if isfield(parameters, 'recording')
    result.recording = data.recording;
end

end



%% ----------------------------------------------------------------------%%

function keepLooping=m_checkConvergence(data, parameters, relativeDifference, time)
keepLooping = true;

%% Number of iteration
if data.iteration == parameters.MaximumIteration
    keepLooping = false;
    if isfield(parameters, 'display')
        fprintf(1, 'Stopping upon reaching maximum iteration (%i).\n', parameters.MaximumIteration);
    end
end

%% Time lapse
if isfield(parameters, 'MaximumTime')
    if time > parameters.MaximumTime
        keepLooping = false;
        if isfield(parameters, 'display')
            fprintf(1, 'Stopping because of the time limit (%fs).\n', parameters.time);
        end
    end
end

%% Relative difference (if allowed to stop the algorithm)
if data.iteration > 1
    if (~isfield(data, 'doNotCheckConvergence')) || (data.doNotCheckConvergence == 0)
        if relativeDifference < parameters.RelativeDifferenceTolerance
            keepLooping = false;
            if isfield(parameters, 'display')
                fprintf(1, 'Stopping because the relative change was too low (%f%%).\n', 100 * relativeDifference);
            end
        end
    end
end

%% key pressed
if isfield(parameters, 'display')
    key = get(gcf,'CurrentCharacter');
    if strcmp(key, 'p') % p: exit key
        keepLooping = false;
        fprintf(1, 'Stopping because key ''%s'' was pressed on the figure.\n', key);
        refresh
        pause(1.5);
    end
end

end


%% ----------------------------------------------------------------------%%

function reference = m_combineFields(reference, addition)
% adds the fields from structure "addition" to structure "reference"
% without replacing preexisting fields of structure reference

if isstruct(addition)
    addFields = listSubfields(addition);
    refFields = listSubfields(reference);
    for k = 1 : length(addFields)
        name = addFields{k};
        if ~sum(ismember(refFields, name))
            reference = setSubfield(reference, name,...
                getSubfield(addition, name));
        end
    end
end

end


%% ----------------------------------------------------------------------%%

function data = m_record(data, parameters)

iter = data.iteration;
if iter == 1 % for initializations
    init = NaN * ones(parameters.MaximumIteration, 1);
end

%check all the fields
recFields = fieldnames(parameters.recording);

for k = 1 : length(recFields)
    val = parameters.recording.(recFields{k})(data);
    if isstruct(val)
        % if a function gave several outputs in a structure
        % save under the outputs name
        valFields = fieldnames(val);
        for j = 1 : length(valFields)
            if iter == 1 % for initialization
                data.recording.(valFields{j}) = init;
            end
            data.recording.(valFields{j})(iter) = val.(valFields{j});
        end
    else
        % otherwise, record under the function name
        if iter == 1 % for initialization
            data.recording.(recFields{k}) = init;
        end
        data.recording.(recFields{k})(iter) = val;
    end
end

end




