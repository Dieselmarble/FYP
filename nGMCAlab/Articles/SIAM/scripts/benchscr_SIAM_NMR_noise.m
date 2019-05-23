clear benchmark
curFolder = pwd;
rootDirectory = 'nGMCAlab';
num = strfind(curFolder, rootDirectory);
cd([curFolder(1:num-1) rootDirectory])
setup_ngmcalab

%% data parameters
benchmark.dataParameters.r = 12;
%number of observations
benchmark.dataParameters.m = 32;
%number of samples
%benchmark.dataParameters.n = 1024;
%shape of A
benchmark.dataParameters.Aalpha = 2;
%activation rate in A
benchmark.dataParameters.Abernoulli = 1;
%noise level
benchmark.dataParameters.dB = linspace(10, 25, 6);
%peaks width
benchmark.dataParameters.width = 4;


%% algorithm parameters

%stopping criteria
benchmark.algorithmParameters.MaximumIteration = 300;
benchmark.algorithmParameters.RelativeDifferenceTolerance = 0.000000001;
benchmark.algorithmParameters.phaseRatio = 0.75;
benchmark.algorithmParameters.S.directSparsity = 1;
benchmark.algorithmParameters.S.qmf = {MakeONFilter('Symmlet', 4)};
benchmark.algorithmParameters.S.L = 3;
%benchmark.algorithmParameters.display = @(data) plot(data.S(1 : 5, :)', 'LineWidth', 2);

%% algorithms
k = 0;

k = k + 1;
benchmark.algorithms(k).name = 'analysis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 1;
benchmark.algorithms(k).displayName = 'rew. ana. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'analysis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 0;
benchmark.algorithms(k).displayName = 'ana. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'synthesis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 0;
benchmark.algorithms(k).displayName = 'syn. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'synthesis_wave_nGMCA_alg';
benchmark.algorithms(k).parameters.S.reweightedL1 = 1;
benchmark.algorithms(k).displayName = 'rew. syn. nGMCA (red. wav.)';

k = k + 1;
benchmark.algorithms(k).name = 'ortho_wave_nGMCA_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'syn. nGMCA (ortho.)';

k = k + 1;
benchmark.algorithms(k).name = 'synthesis_conv_nGMCA_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'syn. nGMCA (conv.)';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCA_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'nGMCA';



%% benchmark settings
clear benchmark_settings
%tag for the name
benchmark.settings.tag = 'SIAM_NMR';
%monte-carlo
benchmark.settings.num_MC = 36;


benchmark.settings.criteriaNames = {'SDR_S', 'SDR_A', 'SIR_S', 'SIR_A', ...
    'SNR_S', 'SNR_A', 'SAR_S', 'SAR_A', 'identificationRate'};
benchmark.settings.drawOptions.criteriaNames = {'SDR_S'};

%eval function
benchmark.settings.criteriaFunction = @(result, reference) NMFevaluation(result, reference, 0);

%data generator
benchmark.settings.dataGenerator = @createReducedRealisticMixtures;


%% Launch
benchmark = performBenchmark(benchmark);

