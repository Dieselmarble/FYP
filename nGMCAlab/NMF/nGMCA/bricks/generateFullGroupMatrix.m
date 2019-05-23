% generateFullGroupMatrix.m - This file is part of nGMCALab.
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

% y = generateFullGroupMatrix(vals, sizes, groups)
% generates a matrix repeating the values in vals for size "sizes"
%
% Inputs:
% - vals: values to repeat.
% - sizes: sizes of the groups.
% - groups: depleted
%
% Output:
% - y: block matrix with value vals(i,j,k) for (i,j,k)-th block of size "sizes" 

function y = generateFullGroupMatrix(vals, sizes, groups)
    if nargin < 3
        groups = [];
    end

    if isempty(groups)
        
        vsizes = size(vals);
        sizes = [sizes, ones(1, 3 - length(sizes))];
        vsizes = [vsizes, ones(1, 3 - length(vsizes))];
        
        
        ysizes = sizes .* vsizes;
        y = zeros(ysizes);
        for i1 = 1 : vsizes(1)
            for i2 = 1 : vsizes(2)
                for i3 = 1 : vsizes(3)
                    y((i1 - 1) * sizes(1) + (1 : sizes(1)),...
                        (i2 - 1) * sizes(2) + (1 : sizes(2)),...
                        (i3 - 1) * sizes(3) + (1 : sizes(3))) = ...
                        vals(i1, i2, i3);
                end
            end
        end
        
    else
        
        y = zeros(size(vals, 1), size(groups, 2), size(groups, 3));
        for r = 1 : size(vals, 1)
            for s = 1 : size(vals, 2)
                y(r, :, :) = y(r, :, :) + vals(r, s) * groups(s, :, :);
            end
        end
        
    end

end

