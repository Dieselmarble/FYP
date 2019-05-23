% mexingRoutine.m - This file is part of the redWaveToolbox.
% This software aims at performing redundant wavelet transformations.
% Copyright 2014 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
% Created on 14/7/2014, last modified on 16/7/2014
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

% mexingRoutine(cpp_name)
% mexes files in the sources directory and move them to the bin directory

function mexingRoutine(cpp_name)

% if cpp_name does not end with '.cpp', add it
if ~strcmp(cpp_name(end - 3 : end), '.cpp')
    cpp_name = [cpp_name, '.cpp'];
end


% mex
mex(['sources/' cpp_name]);

% find mex name (depends on the platform)
name = dir([cpp_name(1 : end - 4) '.mex*']);

% delete previous version if it exists
% (prevents error when the mex was loaded)
if exist(['bin/' name.name],'file') == 3
    delete(['bin/' name.name])
end

% move new version to bin folder
movefile(name.name, 'bin')



end



