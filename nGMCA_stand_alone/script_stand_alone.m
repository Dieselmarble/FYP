%Example of use:
m = 200; %number of observations
n = 200; %number of source samples
r = 5; %number of sources
ref_A = rand(m, r); %non-negative mixtures
ref_S = rand(r ,m) .* (rand(r, m) > 0.9); %non-negative and sparse sources
Y = ref_A * ref_S + 0.1 * randn(m, n); %simulated data
[A, S] = nGMCA(Y, r); %perform decomposition

%same but without text
parameters.verbose = 0;
[A,S] = nGMCA(Y, r, parameters);
%more options are available, check help for more information

%plot
subplot(2, 1, 1); plot(ref_S' / max(max(ref_S)), 'LineWidth', 2);
set(gca, 'FontName', 'TimesNewsRoman');
title('Reference sources', 'FontName', 'TimesNewsRoman');
subplot(2,1,2); plot(S' / max(max(S)), 'LineWidth', 2);
title('Estimated sources (up to scale and permutation invariances)', 'FontName', 'TimesNewsRoman');
set(gca, 'FontName', 'TimesNewsRoman');
% colors do not match because the sources are not reordered.