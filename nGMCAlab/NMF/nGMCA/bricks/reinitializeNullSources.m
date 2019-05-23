% reinitializeNullSources.m - This file is part of nGMCALab.
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
% data = reinitializeNullSources(data,Y,verbose)
% Straightforward reinitialization of a line of S and a column of A
% by picking one column in the residue.
%
% Inputs:
% - data: structure with fields A and S.
% - Y: the data matrix..
% - verbose: 1 for printing information, 0 otherwise.
%
% Output:
% - data : same structure as input, after reinitialization of null sources
%   Non-null sources are left untouched.
function data = reinitializeNullSources(data, Y, verbose)


if nargin < 3
    verbose = 1;
end

if sum(sum(data.S, 2) == 0) > 0 || sum(sum(isnan(data.S))) > 0
    if verbose
        fprintf(1, 'Reinitialization of %i null line(s) in S.\n',...
            sum(sum(data.S, 2) == 0));
    end
    indices = ((sum(data.S, 2) == 0) + sum(isnan(data.S), 2) > 0) > 0;
    [data.A(:, indices), data.S(indices, :)] =...
        m_fastExtractNMF(Y - data.A * data.S, sum(indices));
end

end



function [A, S] = m_fastExtractNMF(residual, r)

if r > 0
    [m, n] = size(residual);
    A = zeros(m, r);
    S = zeros(r, n);
    for i = 1 : r
        residual = residual .* (residual > 0);
        if sum(residual(:)) ~= 0
            % compute square norm of residual to select maximum one
            res2 = sum(residual .* residual);
            j = find(res2 == max(res2));
            if ~isempty(j)
                j = j(1);
                if res2(j) > 0
                    %normalize maximum residual
                    A(:, i) = residual(:, j) / sqrt(res2(j));
                    %compute scalar product with the rest od the residual and keep only
                    %positive coefficients
                    S(i, :) = A(:, i)' * residual;
                    S(i, :) = max(S(i, :), 0);
                    %compute new residual
                    residual = residual - A(:, i) * S(i, :);
                end
            end
        end
    end
else
    A = [];
    S = [];
end
end

