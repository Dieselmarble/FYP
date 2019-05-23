% This script gives an example on how to use on simulated mixtures of
% 2D images the algorithms developped in:
% Jeremy Rapin, Jerome Bobin, Anthony Larue, Jean-Luc Starck. 
% Sparse and Non-negative BSS for Noisy Data,
% IEEE Transactions on Signal Processing, vol. 61, issue 22, p. 5620-5632, 2013. 

setup_ngmcalab; % initialization of the toolbox

%% data settings
clear data_settings parameters
data_settings.r = 10;%number of sources
data_settings.m = 200;%number of observations
data_settings.n = 200;%number of samples
data_settings.Aalpha = 2;%shape of A
data_settings.Abernoulli = 1;%activation rate in A
data_settings.Salpha = 1;%shape of S
data_settings.Sbernoulli = 0.08;%activation rate in S %0.08
data_settings.dB = 10;

%% create the data
gen_seed = 32;%set a seed for repetability with random generators
randGenSeed(gen_seed);
reference = createSparseData(data_settings);
reference = NMFnormalization(reference,'A');

%print exact dB
if isfinite(data_settings.dB)
    fprintf(1, '\n# Noisy data #\n');
    fprintf(1, '%f dB.\n\n', 20 * log10(norm(al(reference.Y), 'fro') / norm(al(reference.N), 'fro')));
end





%% parametrize the algorithms

clear parameters
parameters.Y = reference.Y + reference.N; %input noisy data
parameters.rank = data_settings.r; %input number of sources to recover
comp_seed = 58;%set a seed for repetability with random generators
parameters.reference = reference;%input reference (for oracle algorithms)

%save computation information throughout the iterations
parameters.recording.crit = @(data) NMFevaluation(data, reference);

%display figures during the information (choose at most one)
%parameters.display = @(data) scatter(data.S(1, :), data.S(2, :)); %2D scatter plot
parameters.display = @(data) plot(data.S(1 : 5, :)'); %plot of the sources
%parameters.display = @(data) SDRPlot(data); %reconstruction evaluation

%algorithm constants
parameters.S.tau_MAD = 1;
parameters.MaximumIteration = 500;


%% launch algorithm

% choose an algorithm
%algo = nGMCA_oracle_alg();
%algo = nGMCAhard_alg();
%algo = nGMCAnaive_alg();
algo = nGMCA_alg();
%algo = how_to_plug_in_alg(); %shows how to plug in an external alg in this setting

close;randGenSeed(comp_seed); result = ApplyAlgorithm(algo, parameters); %launch
[criteria, result] = NMFevaluation(result, reference, 1); %evaluate
% after a call to NMFevaluation, sources in "result" are reordered
% according the sources in reference.


% equivalent to nGMCA_alg, one can use the nGMCA stand alone function
parameters.tau_MAD = parameters.S.tau_MAD; %old way to set the parameters
randGenSeed(comp_seed); [A, S] = nGMCA(parameters.Y, parameters.rank, parameters);
result.A = A;
result.S = S;
[criteria, result] = NMFevaluation(result, reference, 1);


%% plot
subplot(2, 1, 1); plot(reference.S' / max(al(reference.S)), 'LineWidth', 2);
set(gca, 'FontName', 'TimesNewsRoman');
title('Reference sources', 'FontName', 'TimesNewsRoman');
subplot(2, 1, 2); plot(result.S' / max(al(result.S)), 'LineWidth', 2);
title('Estimated sources (up to scale invariance)', 'FontName', 'TimesNewsRoman');
set(gca, 'FontName', 'TimesNewsRoman');
% colors match because the sources have been reordered by NMFevaluation

%%
%% NMR spectra
%%

clear data_settings
data_settings.m = 15;
data_settings.n = 1200;
data_settings.dB = 20;
data_settings.Aalpha = 2;
data_settings.Abernoulli = 1;
reference = createRealisticMixtures(data_settings);

% plot one spectrum
num = 4;
plot(reference.ppm,reference.S(num,:)', 'color', 'black', 'LineWidth', 2);
xlabel('ppm')
title(reference.names{num}, 'FontName', 'Times New Roman', 'FontSize', 13);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13, 'YTick', 0, 'Box', 'off', 'TickDir', 'out');

%% perform NMF
clear parameters
parameters.Y = reference.Y + reference.N;
parameters.rank = 15;
parameters.display = @(data) showSignalsDecomposition(data, reference);
parameters.S.tau_MAD = 1;
parameters.MaximumIteration = 500;
parameters.recording.lambda = @(data) data.lambda;

%apply
algo = nGMCAhard_alg();
algo = nGMCA_alg();
close;randGenSeed(comp_seed); result = ApplyAlgorithm(algo,parameters);
[criteria, result] = NMFevaluation(result, reference, 1);


% evolution of lambda throughout the iterations:
close
plot(result.recording.lambda, 'color', 'black', 'LineWidth', 2)
title('Evolution of \lambda throughout the iterations', 'FontName', 'Times New Roman', 'FontSize', 13);
xlabel('iteration', 'FontName','Times New Roman', 'FontSize', 13)
ylabel('\lambda', 'FontName','Times New Roman', 'FontSize',13)
