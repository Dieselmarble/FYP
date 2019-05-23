% dimNorm.m - This file is part of nGMCALab.
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
% n = dimNorm(x,dim,type)
% Provides the norm 1 or 2 or Inf along any given dimension of a matrix
%
% Inputs: 
% - x: the matrix.
% - dim: dimension of the computation.
% - type: type of the norm among 1, 2 or Inf
%
function n = dimNorm(x, dim,norm_type)
    if nargin < 2
        dim = 1;
    end
    if nargin < 3
        norm_type = 2;
    end
    
    if norm_type == 2
        n = sqrt(sum(x .* x, dim));
    else
        if norm_type == 1
            n = sum(abs(x), dim);
        else
            if norm_type == Inf
                n = max(abs(x), [], dim); 
            else
                error('Authorized values for norm_type are: 1, 2 and Inf.\n')
            end
        end
    end
end




