% createSparseData.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13, last modified on 16/7/2014
% 
% This software is governed by the CeCILL  license under French law and
% abiding by the rules of distribution of free software.  You can  use,  
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA,  CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy, 
% modify and redistribute granted by the license,  users are provided only
% with a limited warranty  and the software's author,   the holder of the
% economic rights,   and the successive licensors  have only  limited
% liability. 
%
% In this respect,  the user's attention is drawn to the risks associated
% with loading,   using,   modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software, 
% that may mean  that it is complicated to manipulate,   and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,   more generally,  to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.
%
% 
%
% data = createSparseData(param)
% param with fields: m,  r,  n,  A_alpha,  A_bernoulli,  S_alpha,  S_bernoulli
% generates structure data with field A and S.
% - A of size m x r, absolute value of a bernoulli generalized gaussian
% matrix with activation rate A_bernoulli and shape parameter A_alpha.
% - S of size r x n, absolute value of a bernoulli generalized gaussian
% matrix with activation rate S_bernoulli and shape parameter S_alpha.
% - Y  =  A * S. A is renormalized so that norm(Y) = sqrt(m x n) (variance of
% coefficients of 1).
% - N of size m x n with standard deviation computed in order to obtain the
% specified dB level.
%
% param requires fields:
% - m: number of observations.
% - n: number of samples.
% - r: number of sourcces.
% - Aalpha: shape parameter of the generalized gaussian distribution of A.
% - Abernoulli: activation parameter of the Bernoulli distribution of A.
% - Salpha: shape parameter of the generalized gaussian distribution of S.
% - Sbernoulli: activation parameter of the Bernoulli distribution of S.
% - dB: level of noise ("Inf" for noiseless data).

function data = createSparseData(param)

%create absolute value of the hadamard product of an iid Bernoulli matrix
%and an iid generalized Gaussian matrix
data.A = abs(gengau2(param.Aalpha, param.m, param.r)) .*...
    (rand(param.m, param.r) < param.Abernoulli);
data.S = abs(gengau2(param.Salpha, param.r, param.n)) .*...
    (rand(param.r, param.n) < param.Sbernoulli);

%normalize to have std(Y_ij) = 1
data.Y = data.A * data.S;
coeff = sqrt(param.m * param.n) / norm(al(data.Y), 'fro');
data.A = data.A * coeff;
data.Y = data.Y * coeff;

%compute noise level
if isfinite(param.dB)
    noise = 10.^(-param.dB / 20) / sqrt(param.m * param.n)...
        * norm(al(data.A * data.S), 'fro');
else
    noise = 0;
end

data.N  =  gengau2(2,  size(data.Y, 1),  size(data.Y, 2))  *  noise;

end

