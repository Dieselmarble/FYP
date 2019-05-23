% NMFnormalization.m - This file is part of nGMCALab.
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
% data = NMFnormalization(data,normalized_matrix_str)
% inputs:
% - data with fields A (m x r) and S (r x n).
% - normalized_matrix_str: a string 
% if normalized_matrix_str == 'A': columns of A are normalized and the
% weigths are set on S.
% if normalized_matrix_str == 'S': columns of S are normalized and the
% weigths are set on A.
%
function data = NMFnormalization(data,normalized_matrix_str, norm_type)
% data =normalizeA(data,true)
% inputs:
% - data with fields A (m x r) and S (r x n)
% - true: a scalar 
% if true~=0: columns of A are normalized and the weigths are set on S
% otherwise, normalizes S and sets the scales on A.


if nargin<2
    normalized_matrix_str = 'A';
end
if nargin < 3
    norm_type = 2;
end

%norms
Ns = dimNorm(data.S', 1, norm_type);
Na = dimNorm(data.A, 1, norm_type);
N = Na .* Ns;
act = N > 0;

%normalize active sources and mixtures
if strcmp(normalized_matrix_str, 'S')
    data.S(act, :) = bsxfun(@times, data.S(act, :), 1 ./ Ns(act)');
    data.A(:, act) = bsxfun(@times, data.A(:, act), N(act) ./ Na(act));
else   
    if strcmp(normalized_matrix_str, 'A')
        data.S(act, :) = bsxfun(@times, data.S(act, :), (N(act) ./ Ns(act))');
        data.A(:, act) = bsxfun(@times, data.A(:, act), 1 ./ Na(act));
    else
        error('Field normalized_matrix must be either string ''A'' or ''S''.\n') 
    end
end

%sets to zero inactive sources and mixtures
data.S(~act, :) = 0 * data.S(~act, :);
data.A(:, ~act) = 0 * data.A(:, ~act);

end
    




