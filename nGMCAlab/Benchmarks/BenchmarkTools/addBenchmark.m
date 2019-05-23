% addBenchmark.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
% Created on 19/6/2014, last modified on 16/7/2014
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


% benchmark = addBenchmark(benchmark, additionalBenchmark)
% Combines two benchmarks together
% Criteria must be the same, and no algorithm should be present in both
% benchmarks. No further check is done, the user should make sure that both
% benchmarks have the same settings.
%
% Inputs :
% - benchmark: initial benchmark
% - additionalBenchmark: benchmark to be added

function benchmark = addBenchmark(benchmark, additionalBenchmark)

%% check criteria
c1 = benchmark.settings.criteriaNames;
c2 = additionalBenchmark.settings.criteriaNames;
if ~isequal(c1,c2)
    error('Criteria must be the same.')
end
clear c1 c2;

%%
a1 = {benchmark.algorithms(:).displayName};
a2 = {additionalBenchmark.algorithms(:).displayName};
if length(setdiff(a1, a2)) < length(a1)
    error('One repeated algorithm.')
end
clear a1 a2

%%
benchmark.algorithms = [benchmark.algorithms, additionalBenchmark.algorithms];
%points x n_algo x n_measure x [mean, std, med]]
benchmark.results = cat(2, benchmark.results, additionalBenchmark.results);

end
%ismember

