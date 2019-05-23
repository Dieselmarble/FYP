% performBenchmark.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13, last modified on 16/7/2014
%
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software. You can use,
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http : //www.cecill.info".
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
%
%
% benchmark = performBenchmark(benchmark, displayResults)
% performs a provided benchmark
%
% Inputs :
% - benchmark : structure with fields :
%	- settings : main settings of the algorithm including fields :
%     - tag : to append the name of the created file with a tag (optional)
%     - num_MC : number of Monte-Carlo runs (at least 2)
%     - algorithmNames : cell of call names of the algorithms to test,
%     used with ApplyAlgorithms. A function how_to_plug_in_alg shows how
%     to plugin an external algorithm in these settings.
%     - criteriaFunction : function handle taking as input the output of
%     the algorithm, in order to evaluate it
%     - criteriaNames : cell of names of criteria to evaluate (outputs of
%     the criteriaFunction.
%     - dataGenerator : handle of a function generating the data, taking a
%     structure as input.
%     - drawOptions (optional) : option structure for the drawings using
%     drawBenchmark.
%     - addDataInputs (optional) : functions with prototype parameters =
%     standardDataInputs(parameters, dataParameters, reference)
%     which sets the data into the parameters structure (typically,
%     provides the data Y, the rank r and the reference data if need be)
%	- dataParameters : structure given as input to the dataGenerator. When a
%	(unique) field is defined by a vector instead of a scalar or a 1x1 cell,
%  the benchmarks carries on evaluating the algorithms on each of these
%  values.
%	- algorithmParameters : structure given as input to the algorithms. When a
%	(unique) field is defined by a vector instead of a scalar or a 1x1 cell,
%  the benchmarks carries on evaluating the algorithms on each of these
%  values. For algorithm specific parameters, use the syntax
%  benchmark.algorithmParameters."algorithm name"."parameter name"
% - displayResults : 1 to display results on the fly, 0 otherwise (default :
% 1)
%
% Output : - benchmark : same as input with an additional field "results" of
% size n_points x n_algo x n_criteria x [mean, std, med]]
% This benchmark variable should be provided as input to drawBenchmark
% function.
% A file is saved in the benchmark folder when a benchmarks is completed,
% under the form :
% bench_ < tested variable>_p < number of points>MC < num_MC>_ < date>_ < time>_ < tag>.mat
%
% N.B. :
%  - example of scripts, used in Rapin et al. 2013, are provided in
%   nGMCAlab/Articles/TSP
%  - example of scripts, used in Rapin et al. 2014, are provided in
%   nGMCAlab/Articles/SIAM/scripts
%  - external algorithms can be easily plugged in this frameword as shown
%   in how_to_plug_in_ext_alg.m


function benchmark = performBenchmark(benchmark, displayResults)

% to uniformize different versions :
% define randGenSeed in order to setup random generator
if ~exist('randGenSeed', 'var')
    if ~exist('rng', 'var')
        randGenSeed = @(i) rand('twister', i);
    else
        randGenSeed = @(i) rng(i, 'twister');
    end
end

offrand = floor((2^32 - 1) * rand());

%%
if ~isfield(benchmark.settings, 'drawOptions')
    benchmark.settings.drawOptions = [];
end

if ~isfield(benchmark.settings, 'addDataInputs')
    benchmark.settings.addDataInputs = @standardDataInputs;
end


if nargin < 2
    displayResults = 1;
end


%% get field names and lengths for the dataParameters
dataParameters_names = listSubfields(benchmark.dataParameters);
dataParameters_lengths = zeros(1, length(dataParameters_names));
for k = 1 : length(dataParameters_names)
    dataParameters_lengths(k) = length(getSubfield(benchmark.dataParameters, dataParameters_names{k}));
end

%% get field names and length for algorithmParameters
%create (empty) specific settings for all algorithms
num_algorithms = length(benchmark.algorithms);
algorithmNames = {benchmark.algorithms(:).name};
for k = 1 : num_algorithms
    name = algorithmNames{k};
    if ~exist(name, 'file') %check existence
        error('ERROR : unknown algorithm "%s".', name);
    end
end

% check that there is no two identical display names
algorithmDisplayNames = {benchmark.algorithms(:).displayName};
if length(algorithmDisplayNames) ~= length(unique(algorithmDisplayNames))
    error('ERROR : at least one algorithm is repeated more than once.');
end
clear algorithmDisplayNames;



%% algorithm parameters %does not work with subfields...
algorithmParameters_names = listSubfields(benchmark.algorithmParameters);
algorithmParameters_lengths = zeros(1, length(algorithmParameters_names));
for k = 1 : length(algorithmParameters_names)
    algorithmParameters_lengths(k) = length(getSubfield(benchmark.algorithmParameters, algorithmParameters_names{k}));
end


%% check the varying parameter and create the output name

param_list = [dataParameters_names; algorithmParameters_names];
param_num = [dataParameters_lengths, algorithmParameters_lengths];
numberOfPoints = max(param_num);


if sum(param_num ~= 1) ~= 1 %check that there is only one variable to benchmark
    error('ERROR : exactly one variable should have multiple values.');
else
    
    %tagging
    if isfield(benchmark.settings, 'tag')
        tag = ['_' benchmark.settings.tag];
    else
        tag = [];
    end
    
    i = (param_num ~= 1);
    benchmark.settings.startTime = datestr(now);
    benchVariable = param_list{i};
    file_name = ['bench_' benchVariable '_p' num2str(numberOfPoints) 'MC' num2str(benchmark.settings.num_MC) '_' benchmark.settings.startTime(1 : end-3) tag '.mat'];
    file_name = regexprep(file_name, '-', '');
    file_name = regexprep(regexprep(file_name, ':', ''), ' ', '_');
    fprintf(1, ['Starting evaluation of ' benchVariable ' parameter.\n']);
    fprintf(1, ['Output file will be ' file_name '.\n']);
end

benchmark.settings.benchVariable = benchVariable;
benchmark.settings.numberOfPoints = numberOfPoints;

%% prepare the data

num_criteria = length(benchmark.settings.criteriaNames);


%MC x n_algo x n_measure
MC_results = zeros(benchmark.settings.num_MC, num_algorithms, num_criteria);
%points x n_algo x n_measure x [mean, std, med]]
benchmark.results = NaN * ones(numberOfPoints, num_algorithms, num_criteria, 3);
benchmark.settings.resultsStructure = 'points x n_algo x n_measure x [mean, std, med]';

if benchmark.settings.num_MC < 2
    error('ERROR : you need to set more than 1 Monte-Carlo sampling.');
end

% check if drawing works correctly
if displayResults
    drawBenchmark(benchmark, benchmark.settings.drawOptions);
end


%% launch the benchmark
%set timer to have an average duration time
timer = LoopTimer(numberOfPoints * benchmark.settings.num_MC);
%use a random permutation in order to mix the settings (for a better
%average time approximation)
perm = randperm(numberOfPoints);

for i = 1 : numberOfPoints
    % current point
    p = perm(i);
    
    %% create the data settings 
    dataParameters = []; % fill it one by one, since cells must be removed
    for k = 1 : length(dataParameters_names)
        name = dataParameters_names{k};
        x = getSubfield(benchmark.dataParameters, name);
        if iscell(x)
            dataParameters = setSubfield(dataParameters, name, x{mod(p-1,length(x))+1});
        else
            dataParameters = setSubfield(dataParameters, name, x(mod(p-1,length(x))+1));
        end
    end
    

    
    %% set parameters 
    parameters = []; % fill it one by one, since cells must be removed
    for k = 1 : length(algorithmParameters_names)
        name = algorithmParameters_names{k};
        x = getSubfield(benchmark.algorithmParameters, name);
        if isa(x,'function_handle')
            parameters = setSubfield(parameters, name, x);
        else
            if iscell(x)
                parameters = setSubfield(parameters, name, x{mod(p-1,length(x))+1});
            else
                parameters = setSubfield(parameters, name, x(mod(p-1,length(x))+1));
            end
        end
    end
    
    
    
    %% launch the Monte-Carlo sampling
    for MC = 1 : benchmark.settings.num_MC
        timer.notifyNewIteration()
        fprintf(1, 'Monte-Carlo %i/%i of point %i/%i. (%s remaining)', MC, benchmark.settings.num_MC, i, numberOfPoints, timer.stringRemainingTime());
        
        %make it random
        seed = offrand + MC + 1000 * i;
        randGenSeed(seed);
        
        %create the data
        reference = benchmark.settings.dataGenerator(dataParameters);
        
        % input data into the parameters
        %parameters.Y = reference.Y + reference.N;
        %parameters.reference = reference;
        %parameters.rank = dataParameters.r;
        parameters = benchmark.settings.addDataInputs(parameters, dataParameters, reference);
        
        %% Optimization
        for alg = 1 : num_algorithms
            
            % perform with algorithm :
            name = algorithmNames{alg};
            algo = eval([name '(benchmark.algorithms(' num2str(alg) ').parameters)']);
            if alg > 1
                fprintf(1, [', ' benchmark.algorithms(alg).displayName]);
            else
                fprintf(1, [' ' benchmark.algorithms(alg).displayName]);
            end
            
            % apply
            seed = 1000000 + 1000 * i + MC;
            randGenSeed(seed); %make it random
            result = ApplyAlgorithm(algo, parameters);
            
            % save the results
            criteria = benchmark.settings.criteriaFunction(result, reference);
            for crit = 1 : num_criteria
                MC_results(MC, alg, crit) = criteria.(benchmark.settings.criteriaNames{crit});
            end
            
            % refreshing (to avoid stuck plots)
            refresh;
            pause(0.003);
        end
        fprintf(1, '\n');
    end
    %MC_results
    %compute and store the mean
    benchmark.results(p, : , : , 1) = mean(MC_results);
    %compute and store the standard deviation
    benchmark.results(p, : , : , 2) = std(MC_results);
    %compute the median
    benchmark.results(p, : , : , 3) = median(MC_results);
    
    if displayResults
        drawBenchmark(benchmark, benchmark.settings.drawOptions);
    end
    pause(0.003);
    
    
end
benchmark.settings.endTime = datestr(now);

%% save
if exist('Benchmarks', 'dir')
    save(['Benchmarks/' file_name], 'benchmark');
else
    save(file_name, 'benchmark');
end

end



function parameters = standardDataInputs(parameters, dataParameters, reference)
        parameters.Y = reference.Y + reference.N;
        parameters.reference = reference;
        parameters.rank = dataParameters.r;
end

