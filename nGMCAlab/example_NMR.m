% This script gives an example on how to use on simulated NMR data
% the algorithms developped in:
% Jeremy Rapin, Jerome Bobin, Anthony Larue, Jean-Luc Starck. 
% NMF with Sparse Regularizations in Transformed Domains, 
% SIAM Journal on Imaging Sciences (accepted), 2014.

setup_ngmcalab; % initialization of the toolbox

%% data creation
dataVariables.r = 12;
% number of observations
dataVariables.m = 32;
% shape of A
dataVariables.Aalpha = 2;
%activation rate in A
dataVariables.Abernoulli = 1;
% noise
dataVariables.dB = 20;
%width of the peaks
dataVariables.width = 4;

%create the data
genSeed = 44;
randGenSeed(genSeed);
reference = createReducedRealisticMixtures(dataVariables);
Y = reference.Y + reference.N; %noisy data
imagesc(reference.Y)

%% algorithm parameters
clear parameters
parameters.MaximumIteration = 300;
parameters.rank = dataVariables.r;
parameters.Y = Y;
parameters.verbose = 1;
parameters.phaseRatio = 0.75;

parameters.S.qmf = MakeONFilter('Symmlet', 4);
parameters.S.L = 3;
parameters.S.directSparsity = 1;
parameters.S.LaplacianFilterWidth = 4; % width of the Laplacian filter (for analysis_conv_nGMCA_alg)
%parameters.S.reweightedL1 = 1; % for reweighted L1


% display during computation
parameters.display = @(data) plot(data.S(1 : 5, :)', 'LineWidth', 2); % display source 1 to 5
% compute SDR, SIR, SNR, SAR at each iteration
parameters.recording.crit = @(data) NMFevaluation(data, reference);


%% computation
%choose one algorithm
algo = nGMCA_alg();
algo = analysis_wave_nGMCA_alg();
algo = synthesis_wave_nGMCA_alg();
algo = ortho_wave_nGMCA_alg();
%algo = synthesis_conv_nGMCA_alg();

%launch computation
estSeed = 12; randGenSeed(estSeed); result = ApplyAlgorithm(algo,parameters);


%% analysis of the results
% compute SDR, SIR, SNR, SAR
[criteria, result] = NMFevaluation(result, reference, 1);

% show evolution of the criteria during iterations
% (recorded using parameters.recording.crit = @(data) NMFevaluation(data, reference);)
SDRPlot(result)

% decomposition of the first 4 signals into target, interferences, noise and
% artifacts
showSignalsDecomposition(result, reference)

