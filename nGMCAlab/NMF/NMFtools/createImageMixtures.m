% createImageMixtures.m - This file is part of nGMCALab.
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
% function data = createImageMixtures(param, set)
% Creates simulated mixtures of r = 4 images
%
% Inputs: 
% - param structure with fields
%   - m: number of observations
%   - dB: level of noise ("Inf" for noiseless data)
% - set: set of images to use
%   set 1 is: Lena, peppers, a boat, Barbara (128 x 128, default).
%   set 2 is: bridge, cameraman, parrot, house (256 x 256).

%Output data structure with fields
% - A of size m x 4, the absolute value of a Gaussian matrix.
% - S of size 4 x n for which is line is an image.
% - Y = A * S. A is renormalized so that norm(Y) = sqrt(m x n) (variance of
% coefficients of 1).
% - N of size m x n with standard deviation computed in order to obtain the
% specified dB level.
% - sourceShape: shape of the images ([128, 128] for set 1, and [256, 256]
% for set 2)
function data = createImageMixtures(param, set)

if nargin < 2
    set = 1;
end

if set == 1
    load images4x128x128
else
    load images4x256x256_2
end
numSources = size(Sim, 1);
shape2line = @(S) reshape(S,[size(S,1),numel(S)/size(S,1)]);
%showImages(Sim, sourcesShape)



%% data
data.S = shape2line(Sim);
data.A = abs(rand(param.m, numSources)) .* (rand(param.m, numSources) > 0);
data.Y = data.A * data.S;
data.A = data.A / max(al(data.Y));
data.Y = data.Y / max(al(data.Y));
data.sourcesShape = [size(Sim, 2), size(Sim, 3)];

n = size(data.S, 2);
% compute noise level
if isfinite(param.dB)
    noise = 10.^(-param.dB / 20)...
        / sqrt(param.m * n) * norm(al(data.A * data.S), 'fro');
else
    noise = 0;
end
data.N = randn(size(data.Y)) * noise;

end



