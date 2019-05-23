% NMFevaluation.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13, last modified on 16/7/2014
% 
% This software is governed by the CeCILL  license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
%
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
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
% [criteria, estimate] = NMFevaluation(reference, estimate, display)
% inputs:
% - reference: structure with fields A and S.
% - estimate: structure with fields A and S.
% - display: 1 to display the results, 0 otherwise (default)
%
% outputs:
% - criteria: structure containing the evaluation criteria. Each estimated source
% and reference source is paired one to one and SDR, SIR, SNR and SAR criteria 
% from Vincent et al. 2006 are computed for both the spectra and mixtures
% coefficient. The mean values for the spectra are saved in the fields
% SDR_S, SIR_S, SNR_S and SAR_S; and the ones for the mixtures coefficients
% are saved in SDR_A, SIR_A, SNR_A and SAR_A.
% Additionaly, the criterion identificationRate aims at estimating how
% many sources have been correctly identified (based how much the estimated
% sources is corrupted with interferences from other sources).
% - estimate: the estimate structure with fields A and S reordered considering
% the order of the signals in the reference structure. An additional field 
% "decomposition" is provided with the decompositions of the spectra
% according to Vincent et al. 2006 (target, interferences, noise and
% artifacts)

function [criteria, estimate] = NMFevaluation(estimate, reference, display)


if nargin < 3
    display = 0;
end


%% shortcut
Ar = reference.A;
Ae = estimate.A;
Ae(isnan(Ae)) = 0;%make sure there is no NaN value

Sr = reference.S;
Se = estimate.S;
Se(isnan(Se)) = 0;%make sure there is no NaN value

n_r = size(Sr, 1);
n_e = size(Se, 1);
%% SNRs criteria Vincent et al.
%% reordering
% Pairwise association estimate and reference sources
SDR_S = computeRowSDRmatrix(Se, Sr);

if n_r <= n_e % cases when the number of searched sources is not equal to the number of actual sources
    ind = munkres(-SDR_S);
    reordering = [ind, setdiff(1 : n_e, ind)];

    % reorder
    Se = Se(ind, :);
    Ae = Ae(:, ind);
    estimate.S = estimate.S(reordering, :);
    estimate.A = estimate.A(:, reordering);
else
    ind = munkres(-SDR_S');

    % reorder
    reordering = [ind, setdiff(1 : n_r, ind)];
    Sr = Sr(reordering, :);
    Ar = Ar(:, reordering);
    
end
clear ind SDR_S reordering temp indices;

%% computation on S
[crit, decompositionS] = decompositionCriteria(Se, Sr, reference.N);
criteria.SDR_S = crit.SDR;
criteria.medSDR_S = crit.medSDR;
criteria.SIR_S = crit.SIR;
criteria.SNR_S = crit.SNR;
criteria.SAR_S = crit.SAR;
estimate.decomposition = decompositionS;


clear decompositionS;

%% computation on A

[crit] = decompositionCriteria(Ae', Ar', reference.N');
criteria.SDR_A = crit.SDR;
criteria.medSDR_A = crit.medSDR;
criteria.SIR_A = crit.SIR;
criteria.SNR_A = crit.SNR;
criteria.SAR_A = crit.SAR;


%% identification rate
% heuristic aiming at finding the number of sources
% which are identifiable
Sref = diag(1 ./ dimNorm(Sr,2)) * Sr;
Sest = diag(1 ./ dimNorm(Se,2)) * Se;

% comute the coefficients of the projection of the estimate sources on
% the reference sources
M = Sest / Sref;
% projection coeff on the target sources
coeff = diag(M);
% maximum of the projection coeff on the non-target sources
%M = abs(M - [diag(diag(M)), zeros(n_e, max(0, n_r - n_e))]);
M(:, 1 : min(n_e, n_r)) = abs(M(:, 1 : min(n_e, n_r)) - diag(diag(M)));
interf = max(M, [], 2);
% sources are considered as identifiable if the target coeff is
% at least twice as large as the other projection coefficients
criteria.identificationRate = mean(coeff > 2 * interf);


%% display the results if need be
if display
    fprintf(1, 'Results of the reconstruction:\n');
    fprintf(1, 'Decomposition criteria on S:\n');
    fprintf(1, '   - Mean SDR: %f.\n', criteria.SDR_S);
    fprintf(1, '   - Mean SIR: %f.\n', criteria.SIR_S);
    fprintf(1, '   - Mean SNR: %f.\n', criteria.SNR_S);
    fprintf(1, '   - Mean SAR: %f.\n', criteria.SAR_S);
    fprintf(1, '   - Identification rate: %f.\n', criteria.identificationRate);

    if sum(size(Ae)==size(Ar))==2
        fprintf(1,'Decomposition criteria on A:\n');
    fprintf(1, '   - Mean SDR: %f.\n', criteria.SDR_A);
    fprintf(1, '   - Mean SIR: %f.\n', criteria.SIR_A);
    fprintf(1, '   - Mean SNR: %f.\n', criteria.SNR_A);
    fprintf(1, '   - Mean SAR: %f.\n', criteria.SAR_A);
    end

    fprintf(1, '\n');
    fprintf(1, 'Identification rate: %f.\n', criteria.identificationRate);

    fprintf(1, '\n');
end


end


function MSDR = computeRowSDRmatrix(Y, X)
%X is a reference of row signals
%Y is an estimate of row signals
%MSNR_ij is the SDR between the i-th row of X with the j-th row of Y

%normalize the reference
X = diag(1 ./ dimNorm(X,2)) * X;

n_x = size(X, 1);
n_y = size(Y, 1);
MSDR = zeros(n_x, n_y);

for n = 1 : n_x
    targets = (Y * X(n, :)') * X(n, :);
    diff = Y - targets;
    
    norm_diff = max(sum(diff .* diff,2), eps);
    norm_targets = max(eps, sum(targets .* targets,2));
    MSDR(n,:) = -10 * log10(norm_diff ./ norm_targets)';
end

end

%input line signals matrices
function [criteria, decomposition] = decompositionCriteria(Se, Sr, noise)
%Se: estimate
%Sr: reference
r = min(size(Se, 1), size(Sr, 1)); %there may be more reference sources than estimated sources


%normalize reference
Sr = diag(1 ./ dimNorm(Sr,2)) * Sr;

%compute the projections
pS = (Se / Sr) * Sr;
SN = [Sr; noise];
pSN = (Se / SN) * SN;

%targets
t = diag(diag(Se(1 : r, :) * Sr(1 : r, :)')) * Sr(1 : r, :); %Sr is normalized
%interferences
i = pS - t;
%noise
n = pSN - pS;
%artifacts
a = Se - pSN;


%SDR: source to distortion ratio
num = t;
den = i + n + a;
norm_num_2 = sum(num .* num, 2);
norm_den_2 = sum(den .* den, 2);
criteria.SDR = mean(10 * log10(max(norm_num_2, eps) ./ max(norm_den_2, eps)));
criteria.medSDR = median(10 * log10(max(norm_num_2,eps)./max(norm_den_2, eps)));

%SIR: source to interferences ratio
num = t;
den = i;
norm_num_2 = sum(num .* num, 2);
norm_den_2 = sum(den .* den, 2);
criteria.SIR = mean(10 * log10(max(norm_num_2, eps) ./ max(norm_den_2, eps)));

%SNR: source to noise ratio
if norm(al(noise))>0 %only if there is noise
num = t + i;
den = n;
norm_num_2 = sum(num .* num, 2);
norm_den_2 = sum(den .* den, 2);
criteria.SNR = mean(10 * log10(max(norm_num_2, eps)./max(norm_den_2, eps)));
else
    criteria.SNR = NaN;
end

%SAR: sources to artifacts ratio
if (size(noise, 1) + size(Sr, 1)) < size(Sr, 2) %if noise + sources form a base, there is no "artifacts"
    num = t + i + n;
    den = a;
    norm_num_2 = sum(num .* num, 2);
    norm_den_2 = sum(den .* den, 2);
    criteria.SAR = mean(10 * log10(max(norm_num_2, eps) ./ max(norm_den_2, eps)));
else
    criteria.SAR = NaN;
end


%% Decomposition
decomposition.target = t;
decomposition.interferences = i;
decomposition.noise = n;
decomposition.artifacts = a;

end




