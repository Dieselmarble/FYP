clear benchmark
curFolder = pwd;
rootDirectory = 'nGMCAlab';
num = findstr(curFolder, rootDirectory);
cd([curFolder(1 : num - 1) rootDirectory])
setup_ngmcalab


%% data settings

%number of sources
benchmark.dataParameters.r = 5:5:35;
%number of observations
benchmark.dataParameters.m = 200;
%number of samples
benchmark.dataParameters.n = 200;
%shape of A
benchmark.dataParameters.Aalpha = 2;
%activation rate in A
benchmark.dataParameters.Abernoulli = 1;
%shape of S
benchmark.dataParameters.Salpha = 1;
%activation rate in S
benchmark.dataParameters.Sbernoulli = 0.3;
%noise level
benchmark.dataParameters.dB = 15;


%% algorithm parameters

%thresholding coefficient
benchmark.algorithmParameters.S.tau_MAD = 1;

%stopping criteria
benchmark.algorithmParameters.MaximumIteration = 500;
benchmark.algorithmParameters.RelativeDifferenceTolerance = 0.000000001;


%% algorithms
k = 1;
benchmark.algorithms(k).name = 'nGMCA_oracle_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'oracle';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCA_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'nGMCA^S';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCAhard_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'nGMCA^H';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCAnaive_alg';
benchmark.algorithms(k).parameters = [];
benchmark.algorithms(k).displayName = 'nGMCA^{naive}';


%% benchmark settings
%tag for the name
benchmark.settings.tag = 'IEEE15dBa30';
%monte-carlo
benchmark.settings.num_MC = 12;%48
%algorithms
benchmark.settings.algorithmNames = {'nGMCA_oracle_alg', 'nGMCA_alg', 'nGMCAhard_alg', 'nGMCAnaive_alg'};
%criteria
benchmark.settings.criteriaFunction = @(result, reference) NMFevaluation(result, reference);
benchmark.settings.criteriaNames = {'SDR_S', 'SDR_A', 'SIR_S', 'SIR_A', 'SNR_S', 'SNR_A', 'identificationRate'};
%data generator
benchmark.settings.dataGenerator = @createSparseData;
%draw options
benchmark.settings.drawOptions.criteriaNames = {'SDR_S', 'SDR_A'};


%% Launch
benchmark = performBenchmark(benchmark);

