clear benchmark
curFolder = pwd;
rootDirectory = 'nGMCAlab';
num = strfind(curFolder, rootDirectory);
cd([curFolder(1 : num-1) rootDirectory])
setup_ngmcalab

%% data parameters
%number of observations
benchmark.dataParameters.r = 4;
benchmark.dataParameters.m = 32;
%noise level
benchmark.dataParameters.dB = linspace(25, 40, 5);



%% algorithm parameters

%stopping criteria
benchmark.algorithmParameters.MaximumIteration = 200;
benchmark.algorithmParameters.RelativeDifferenceTolerance = 0.000000001;
benchmark.algorithmParameters.phaseRatio = 0.75;
benchmark.algorithmParameters.A.noL2CoarseScale = 1;
benchmark.algorithmParameters.S.sourcesShape = {[128, 128]};
benchmark.algorithmParameters.S.directSparsity = 0;
benchmark.algorithmParameters.S.qmf = {MakeONFilter('Daubechies', 4)};%{MakeONFilter('Symmlet', 4)};
benchmark.algorithmParameters.S.L = 3;

%% algorithms
k = 0;

k = k + 1;
benchmark.algorithms(k).name = 'analysis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 1;
benchmark.algorithms(k).parameters.A.noL2CoarseScale = 1;
benchmark.algorithms(k).displayName = 'rew. ana. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'analysis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 0;
benchmark.algorithms(k).parameters.A.noL2CoarseScale = 1;
benchmark.algorithms(k).displayName = 'ana. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'synthesis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 0;
benchmark.algorithms(k).parameters.A.noL2CoarseScale = 1;
benchmark.algorithms(k).displayName = 'syn. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'synthesis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 1;
benchmark.algorithms(k).parameters.A.noL2CoarseScale = 1;
benchmark.algorithms(k).displayName = 'rew. syn. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCA_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'nGMCA';



%% benchmark settings
clear benchmark_settings
%tag for the name
benchmark.settings.tag = 'SIAM_image';
%monte-carlo
benchmark.settings.num_MC = 24;


benchmark.settings.criteriaNames = {'SDR_S', 'SDR_A', 'SIR_S', 'SIR_A', ...
    'SNR_S', 'SNR_A', 'SAR_S', 'SAR_A', 'identificationRate'};
benchmark.settings.drawOptions.criteriaNames = {'SDR_S'};

%eval function
benchmark.settings.criteriaFunction = @(result, reference) NMFevaluation(result, reference, 0);

%data generator
benchmark.settings.dataGenerator = @createImageMixtures;


%% Launch
benchmark = performBenchmark(benchmark);

