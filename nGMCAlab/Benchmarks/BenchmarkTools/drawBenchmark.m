% drawBenchmark.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
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
% fighandle = drawBenchmark(benchmark, options)
% Function which draws benchmarks created with performBenchmark.m
%
% Inputs :
% - benchmark: benchmark structure to be drawn
% - options: optional structure. Possible fields include :
%  - pointSelection: index of points to be drawn if not all of them
%  - LineWidth: line width of the curves (default : 2)
%  - fontSize: font size of the text (default : 12)
%  - algorithmNames: selection of algorithms to draw if not all of them
%  - criteriaNames: selection of criteria to draw (in subplots) if not all
%  of them.
%  - mode: drawing mode, mean of the values with no error bar (0, default),
%  with error bar (1) or median (2).
%  - colorList: lists of colors for the curves.
%  - markerList: lists of markers for the curves.
%  - lineList: lists of line styles for the curves.

function fighandle = drawBenchmark(benchmark, options)
clf
fontName = 'Times New Roman';
str = computer;
if strcmp(str(1:3),'MAC')
    fontName = 'Helvetica';
end

%%
if nargin < 2
    options = [];
end

if ~isfield(options, 'mode')
    options.mode = 1;
    % 0 : mean, no error bar
    % 1 : mean, with error bar
    % 2 : median
end

if ~isfield(options, 'pointSelection')
    pointSelection = 1 : benchmark.settings.numberOfPoints;
else
    pointSelection = options.pointSelection;
end

if ~isfield(options, 'LineWidth')
    options.LineWidth = 2;
end



if ~isfield(options, 'fontSize')
    options.fontSize = 12;
end

if ~isfield(options, 'coordinates')
    field_names = listSubfields(benchmark.algorithmParameters);
    bench_var = benchmark.settings.benchVariable;
    if sum(ismember(field_names, bench_var))
        coordinates = getSubfield(benchmark.algorithmParameters, bench_var);
    else
        coordinates = getSubfield(benchmark.dataParameters, bench_var);
    end
else
    if isfield(options, 'coordinatesName')
        benchmark.settings.benchVariable = options.coordinatesName;
        coordinates = options.coordinates;
    else
        error('In order to change the coordinate system, \n you need to input a coordinateName in the options structure.\n')
    end
end


if ~isfield(options, 'algorithmNames')
    options.algorithmNames = {benchmark.algorithms(:).displayName};
end



if isfield(options, 'referenceAlgorithm')
    benchmark = setReferenceAlgorithm(options.referenceAlgorithm, benchmark);
end



if ~isfield(options, 'criteriaNames')
    options.criteriaNames = benchmark.settings.criteriaNames;
end

if ~isfield(options, 'colorList')
    colorList = {'blue', 'red', 'green', 'magenta', 'black', 'cyan', [0.4, 0.3, 0.7], [1 0.5 0]};
else
    colorList = options.colorList;
end


if ~isfield(options, 'markerList')
    markerList = {'x', 'o', 'p', '^', 's', 'd', 'v', '>', ' < ', '+', '*', 'h'};
else
    markerList = options.markerList;
end


if ~isfield(options, 'lineList')
    lineList = {'--', '-.'};
else
    lineList = options.lineList;
end


%% draw
num_algorithms = length(options.algorithmNames);
num_criteria = length(options.criteriaNames);


%parameters
mx = min(coordinates);
Mx = max(coordinates);
mx = mx - (Mx - mx) / 25;
Mx = Mx + (max(coordinates) - min(coordinates)) / 25;


for i = 1 : num_criteria
    %find the index of the criteria
    crit = find(ismember(benchmark.settings.criteriaNames, options.criteriaNames{i}) == 1);
    if isempty(crit)
        error('ERROR : criterium %s not found.\n', options.criteriaNames{i})
    end
    
    subplot(num_criteria, 1, i)
    for j = 1 : num_algorithms
        
        %find the index of the algorithm
        alg = find(ismember({benchmark.algorithms(:).displayName}, options.algorithmNames{j}) == 1);
        if isempty(alg)
            error('ERROR : algorithm ''%s'' not found.\n', options.algorithmNames{j})
        end
        
        %prepare the values
        col = colorList{mod(j - 1, numel(colorList)) + 1};
        coord = coordinates(pointSelection);
        switch options.mode
            case 2
                vals = benchmark.results(pointSelection, alg, crit, 3)';
            case 3
                vals = benchmark.results(pointSelection, alg, crit, 2)';
            otherwise
                vals = benchmark.results(pointSelection, alg, crit, 1)';
        end
%         if options.mode == 2
%             vals = benchmark.results(pointSelection, alg, crit, 3)';
%         else
%             vals = benchmark.results(pointSelection, alg, crit, 1)';
%         end
        
        
        
        marker = markerList{mod(j - 1, numel(markerList)) + 1};
        line = lineList{mod(j - 1, numel(lineList)) + 1};
        markersize = 10;
        if marker == 'x'
            markersize = 12;
        end
        
        
        if options.mode == 1
            valstds = benchmark.results(pointSelection, alg, crit, 2)';
            errorbar(coord, vals, valstds, line, 'color', col, 'LineWidth', options.LineWidth, 'Marker', marker, 'MarkerSize', markersize);
        else
            plot(coord, vals, line, 'color', col, 'LineWidth', options.LineWidth, 'Marker', marker, 'MarkerSize', markersize);
        end
        hold on
    end
    v = axis();
    axis([mx, Mx, v(3 : 4)]);
    xlabel(benchmark.settings.benchVariable, 'FontSize', options.fontSize, 'FontName', 'Times New Roman')
    ylabel(benchmark.settings.criteriaNames{crit}, 'FontSize', options.fontSize, 'FontName', 'Times New Roman')
    
    set(gca, ...
        'Box'     , 'off'   , ...
        'TickDir'   , 'out'   , ...%'TickLength' , [1 1]) , ...
        'XMinorTick' , 'on'   , ...
        'YMinorTick' , 'on'   , ...
        'YGrid'    , 'on'   , ...
        'XColor'   , [.15 .15 .15], ...
        'YColor'   , [.15 .15 .15], ...
        'LineWidth'  , 1, ...
        'fontsize', options.fontSize, ...
        'FontName', fontName);
end


legend_text = ['legend(''' options.algorithmNames{1} ''''];
for k = 2 : num_algorithms
    legend_text = [legend_text ', ''' options.algorithmNames{k} ''''];
end
legend_text = [legend_text ', 0);'];
eval(['AX = ', legend_text]);


LEG = findobj(AX, 'type', 'text');
set(LEG, 'FontSize', options.fontSize, 'FontName', fontName);

fighandle = gcf;


end




function benchmark = setReferenceAlgorithm(reference_algorithm, benchmark)
%% Setting an algorithm as reference
%%
% find the index
index = 0;
for k = 1 : length(benchmark.settings.algorithmNames)
    if strcmp(benchmark.settings.algorithmNames{k}, reference_algorithm)
        index = k;
    end
end
if index
    [n_points, n_algo, n_measures, mean_and_std] = size(benchmark.results);
    for i = 1 : n_points
        for k = 1 : n_measures
            ref = benchmark.results(i, index, k, 1);
            for j = 1 : n_algo
                benchmark.results(i, j, k, 1) = benchmark.results(i, j, k, 1)-ref;
            end
        end
    end
end
end

