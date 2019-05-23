% computeWaveScaleFun.m - This file is part of nGMCALab.
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

% Dimension 1 and only dimension 1 must be a non-wavelet dimension

% output = computeWaveScaleFun(size_X, Xw, fun, groups)
% computes a specified function on each scale of Xw (in the wavelet
% domain), and on each source if there are several sources.
% To be used alongside with generateFullGroupMatrix which can create a
% matrix of the size of Xw, which evaluates at the value computed for the
% scale in which the coefficient is.
%
% Inputs:
%  - size_X: size of the matrix X in the direct domain.
%  - Xw: matrix in the wavelet domain.
%  - fun: function to compute on each scale of the transform.
%  - groups: to specify groups instead of each scale of the transform
%  (depleted).
%
% Output:
%  - output: each element is the value of the function for a specified scale
%  and source.


function output = computeWaveScaleFun(size_X, Xw, fun, groups)
if nargin < 4
    groups = [];
end

size_Xw = size(Xw);
size_Xw = [size_Xw, ones(1, 3 - length(size_Xw))];
size_X = [size_X, ones(1, 3 - length(size_X))];
if size_Xw(1) ~= size_X(1)
    error('Dimension 1 and only 1 must be a non-wavelet dimension')
end
if length(size(Xw)) > 3
    error('Only up to 2D sources allowed')
end

num_scales = size_Xw ./ size_X;
n1 = size_X(1);
n2 = size_X(2);
n3 = size_X(3);

if ~isempty(groups)
    output = zeros(n1, size(groups, 1));
    for r = 1 : n1
        xw = Xw(r, :, :);
        for s = 1 : size(groups, 1)
            output(r, s) =  fun(xw(groups(s, :, :)));
        end
    end
else
    output = zeros([n1, num_scales(2 : end)]);
    for r = 1 : n1
        for i2 = 1 : num_scales(2)
            for i3 = 1 : num_scales(3)
                output(r, i2, i3) = fun(Xw(r,...
                    n2 * (i2 - 1) + (1 : n2),...
                    n3 * (i3 - 1) + (1 : n3)));
            end
        end
    end
end
end
