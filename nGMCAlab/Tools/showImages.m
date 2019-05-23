% showImages.m - This file is part of nGMCALab.
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

%
% showImages(S, imSize, indices)
% Display function for multispectral image data.
% Inputs: 
% - S: image data of size [N, imSize(1), imSize(2)] where N is the number
% of images.
% - imSize: size of the 2D images
% - indices: indices of the images to be displayed, in [1,..,N]
%
function showImages(S, imSize, indices)

if nargin < 2
    imSize = [sqrt(size(S, 2)), sqrt(size(S, 2))];
end
if nargin < 3
    indices = 1 : size(S, 1);
end


if size(S, 3) == 1%if in line shape
    Sim = permute(reshape(S, [size(S, 1), imSize]), [2 3 1]);
else
    Sim = permute(S, [2 3 1]);
end

% sizes of the subplots
num = length(indices);
sx = floor(sqrt(num) + 0.99999);
sy = floor(num / sx + 0.99999);


for k = 1 : num
    subplot(sy, sx, k)
    m_plotImage(Sim, indices(k))
end

end





function m_plotImage(Sim, ind)
imagesc(Sim(:, :, ind)); colormap(linspace(0, 1, 256)' * [1, 1, 1]);
set(gca,'DataAspectRatio',[1 1 1])
title(['source #' num2str(ind)], 'FontName', 'Times New Roman')
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'XColor'      , 'black', ...
    'YColor'      , 'black', ...
    'YTick'       , [], ...
    'XTick'       , [], ...
    'LineWidth'   , 1         ,...
    'FontName','Times New Roman');
hold off
end








