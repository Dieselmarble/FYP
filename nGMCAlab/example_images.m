% This script gives an example on how to use on simulated mixtures of
% 2D images the algorithms developped in:
% Jeremy Rapin, Jerome Bobin, Anthony Larue, Jean-Luc Starck. 
% NMF with Sparse Regularizations in Transformed Domains, 
% SIAM Journal on Imaging Sciences (accepted), 2014.

setup_ngmcalab; % initialization of the toolbox

%% create data
dataParameters.m = 16; % number of observations
dataParameters.dB = 30; % noise level
genSeed = 1024;
randGenSeed(genSeed);
reference = createImageMixtures(dataParameters); % create data
sourcesShape = reference.sourcesShape; %2D shape of images
numSources = size(reference.S, 1);
Y = reference.Y + reference.N; % noisy data. Each observation is a row of Y (images are vectorized)
showImages(reference.Y + reference.N, sourcesShape) % show data

%% parametrize
% (all parameters have default value if not manually set, except Y and rank)
clear parameters
parameters.rank = numSources; % number of sources to recover
parameters.Y = Y; % noisy data
parameters.phaseRatio = 0.65; %ratio of iterations before refinement phase
parameters.RelativeDifferenceTolerance = 0.000000001;
parameters.MaximumIteration = 200;
parameters.verbose = 1;

% parameters for the update of A
parameters.A.MaximumIteration = 120; % subiterations for A
parameters.A.RelativeDifferenceTolerance = 0.0000000001;
parameters.A.noL2CoarseScale = 1;

% parameters for the update of S
parameters.S.tau_MAD = 2; %thresholding at tau_MAD * sigma
parameters.S.verbose = 1;
parameters.S.MaximumIteration = 24; % subiterations for S
parameters.S.L = 3; % number of wavelet scales
parameters.S.qmf = MakeONFilter('Symmlet', 4); % wavelet filter
parameters.S.sourcesShape = sourcesShape; % 2D shape of the images
parameters.S.directSparsity = 0; % whether the images are sparse in the direct domain or not
parameters.S.reweightedL1 = 1; % use reweigthed L1

% function to display at each iteration
parameters.display = @(data) showImages(data.S, sourcesShape); % show current images
% function(s) to compute at each iteration
parameters.recording.crit = @(data) NMFevaluation(data, reference); % compute criteria at each iteration
%%
%algo = synthesis_wave_nGMCA_alg();
algo = analysis_wave_nGMCA_alg();

% computation
estSeed = 12; close;randGenSeed(estSeed); result = ApplyAlgorithm(algo, parameters);%seed=10;

%% evaluation of the results

% compute SDR, SIR, SNR, SAR
[criteria, result] = NMFevaluation(result,reference,1);

% show result images
showImages(result.S)

% show SDR, SIR, SNR, SAR during all iterations
% (as computed by parameters.recording.crit)
SDRPlot(result)

% show image decomposition into interferences, noise and artifacts
showImageDecomposition(result, sourcesShape)
