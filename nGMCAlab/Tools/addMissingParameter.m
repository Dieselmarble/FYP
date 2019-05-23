% addMissingParameter.m - This file is part of nGMCALab.
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
% parameters = addMissingParameter(parameters, field_name, default_value, function_name)
% Inputs:
% - parameters: a structure.
% - field_name: string name of the field to add. Subfields can be set by
% using a dot in the name.
% - default_value: value of the field to add.
% - function_name: name of the calling algorithm (for verbose mode).
%
% Output:
% parameters: parameters structure with field field_name to default_value
% if it was missing, unchanged otherwise.
% If parameters.verbose ~= 0, displays notifications when adding a missing field.
%

function parameters = addMissingParameter(parameters, field_name, default_value, function_name, verbose)

%list all the fields
allFields = listSubfields(parameters);

if nargin < 3
    if ~sum(ismember(allFields, field_name))
        error('Parameter "%s" is required.\n', field_name);
    end
else
    if nargin < 4
        function_name = 'unknown function'; 
    end
    if nargin < 5
        verbose = 0; 
        if isfield(parameters, 'verbose')
            verbose = parameters.verbose; 
        end
    end

    if ~sum(ismember(allFields, field_name))
        parameters = setSubfield(parameters, field_name, default_value);
        if verbose > 0
            if isscalar(default_value) && ~isstruct(default_value)
                fprintf(1,'In %s, parameter "%s" automatically set to %f.\n',...
                    function_name, field_name, default_value);
            else
                fprintf(1,'In %s, parameter "%s" automatically set to default value.\n',...
                    function_name, field_name);
            end
        end
    end
end

end



