clear benchmark
curFolder = pwd;
rootDirectory = 'nGMCAlab';
num = findstr(curFolder, rootDirectory);
cd([curFolder(1 : num - 1) rootDirectory])
setup_ngmcalab

%% data settings

benchmark.dataParameters.r = 15;
%number of observations
benchmark.dataParameters.m = 15;
%number of samples
benchmark.dataParameters.n = 1200;
%shape of A
benchmark.dataParameters.Aalpha = 2;
%activation rate in A
benchmark.dataParameters.Abernoulli = 1;
%noise level
benchmark.dataParameters.dB = 10 : 5 : 30;


%% algorithm parameters

%thresholding coefficient
benchmark.algorithmParameters.tau_MAD = 2;

%stopping criteria
benchmark.algorithmParameters.MaximumIteration = 500;
benchmark.algorithmParameters.RelativeDifferenceTolerance = 0.000000001;


%% algorithms
k = 1;
benchmark.algorithms(k).name = 'nGMCA_oracle_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 2;
benchmark.algorithms(k).displayName = 'oracle';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCA_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 2;
benchmark.algorithms(k).displayName = 'nGMCA^S';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCAhard_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 2;
benchmark.algorithms(k).displayName = 'nGMCA^H';

k = k + 1;
benchmark.algorithms(k).name = 'nGMCAnaive_alg';
benchmark.algorithms(k).parameters.S.tau_MAD = 2;
benchmark.algorithms(k).displayName = 'nGMCA^{naive}';


%% benchmark settings

%tag for the name
benchmark.settings.tag = 'IEEEappli';
%monte-carlo
benchmark.settings.num_MC = 36;
%criteria
benchmark.settings.criteriaFunction = @(result, reference) NMFevaluation(result, reference);
benchmark.settings.criteriaNames = {'SDR_S', 'SDR_A', 'SIR_S', 'SIR_A', 'SNR_S', 'SNR_A', 'identificationRate'};
benchmark.settings.drawOptions.criteriaNames = {'SDR_S', 'SDR_A'};

%data generator
benchmark.settings.dataGenerator = @createRealisticMixtures;


%% Launch
benchmark = performBenchmark(benchmark);

