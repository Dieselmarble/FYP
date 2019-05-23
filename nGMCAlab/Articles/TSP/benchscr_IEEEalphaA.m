clear benchmark
curFolder = pwd;
rootDirectory = 'nGMCAlab';
num = findstr(curFolder, rootDirectory);
cd([curFolder(1 : num - 1) rootDirectory])
setup_ngmcalab


%% data settings

%number of sources
benchmark.dataParameters.r = 35;
%number of observations
benchmark.dataParameters.m = 200;
%number of samples
benchmark.dataParameters.n = 200;
%shape of A
benchmark.dataParameters.Aalpha = [0.5000 0.7000 0.9000 1.2000 1.6000 2.0000 3.0000 4.0000];
%activation rate in A
benchmark.dataParameters.Abernoulli = 1;
%shape of S
benchmark.dataParameters.Salpha = 1;
%activation rate in S
benchmark.dataParameters.Sbernoulli = 0.8;
%noise level
benchmark.dataParameters.dB = Inf;


%% algorithm parameters

%thresholding coefficient
benchmark.algorithmParameters.S.tau_MAD = 0;

%stopping criteria
benchmark.algorithmParameters.MaximumIteration = 500;
benchmark.algorithmParameters.RelativeDifferenceTolerance = 0.000000001;


%% algorithms
k = 1;
benchmark.algorithms(k).name = 'nGMCA_oracle_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 0;
benchmark.algorithms(k).displayName = 'oracle';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCA_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 0;
benchmark.algorithms(k).displayName = 'nGMCA^S';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCAhard_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 0;
benchmark.algorithms(k).displayName = 'nGMCA^H';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCAnaive_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 0;
benchmark.algorithms(k).displayName = 'nGMCA^{naive}';


%% benchmark settings

%tag for the name
benchmark.settings.tag = 'IEEEr35a80';
%monte-carlo
benchmark.settings.num_MC = 48;
%criteria
benchmark.settings.criteriaFunction = @(result, reference) NMFevaluation(result, reference);
benchmark.settings.criteriaNames = {'SDR_S', 'SDR_A', 'SIR_S', 'SIR_A', 'SNR_S', 'SNR_A', 'identificationRate'};
benchmark.settings.drawOptions.criteriaNames = {'SDR_S', 'SDR_A'};

%data generator
benchmark.settings.dataGenerator = @createSparseData;


%% Launch
benchmark = performBenchmark(benchmark);

