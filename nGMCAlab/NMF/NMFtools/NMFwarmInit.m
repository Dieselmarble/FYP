% NMFwarmInit.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : J�r�my Rapin (jeremy.rapin@cea.fr)
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
% function data = NMFwarmInit(Y, r)
% Initialization of the data structure with fields A and S.
% The initialization consists in alternating between a couple of unconstrained
% (fast) updates and constrained (more precise) updates.
%
% Inputs:
% - Y: data matrix.
% - r: number of sources.
% Output: data with fields A and S such that Y \approx A * S
%
function data = NMFwarmInit(Y, r)


    m = size(Y, 1);
    data.A = rand(m, r);
    options.MaximumIteration = 50;
    options.RelativeDifferenceTolerance = 0;
    
    for k=1:2
        
        %% unconstrained updates
        i = sum(data.A) > 0;
        data.S(i, :) = max(data.A(:, i) \ Y, 0);
        i = sum(data.S,2)>0;
        data.A(:, i) = max(Y / data.S(i, :), 0);


        %% constrained updates
        data.S = nonnegativeSparseUpdate(data.S, data.A' * data.A,...
            data.A' * Y, 0, options);
        % reinitialize S if need be
        data = reinitializeNullSources(data, Y, 0);
        data.A = nonnegativeSparseUpdate(data.A', data.S * data.S',...
        	data.S * Y', 0, options)';
    end
end

