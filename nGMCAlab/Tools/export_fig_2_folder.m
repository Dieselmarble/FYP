% export_fig_2_folder.m - This file is part of nGMCALab.
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


% export_fig_2_folder(fig_name, folder_name)
% Figure exportation tool, using the export_fig toolbox. Convenient for
% publication since it saves the file in a publishable eps format in the
% folder of your choice, and deletes temporary pdf files created by
% pdflatex.
%
% Inputs:
% - fig_name: name to be given to the eps file of the figure.
% - folder_name: name of the folder in which to put the file.

function export_fig_2_folder(fig_name, folder_name)

if ~exist('export_fig', 'file')
    fprintf(1, 'export_fig_2_folder requires the export_fig toolbox. Skipping the instruction.\n')
else
    %export to eps only
    if ~strcmp(fig_name(end - 3 : end), '.eps')
        fig_name = [fig_name, '.eps'];
    end
    if ~strcmp(folder_name(end), '\')
        if ~strcmp(folder_name(end), '/')
            folder_name = [folder_name, '\'];
        end
    end
    
    %check for correct disposition
    fprintf(1, 'Check correct disposition and press a key.\n');
    pause()
    
    
    
    %remove pdf files
    if exist([folder_name,fig_name(1 : end - 3), 'pdf'], 'file')
        delete([folder_name, fig_name(1 : end - 3), 'pdf']);
    end
    if exist([folder_name, fig_name(1 : end - 4), '-eps-converted-to.pdf'], 'file')
        delete([folder_name, fig_name(1 : end - 4), '-eps-converted-to.pdf'])
    end
    %and previous eps
    if exist([folder_name, fig_name], 'file')
        delete([folder_name, fig_name])
    end
    
    %export
    export_fig temp.eps -transparent
    movefile('temp.eps', fig_name);
    movefile(fig_name, folder_name);
    
end

end




