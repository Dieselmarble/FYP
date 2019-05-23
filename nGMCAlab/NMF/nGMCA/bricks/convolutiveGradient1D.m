% convolutiveGradient1D.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
% Created on 7/7/2014, last modified on 16/7/2014
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software. You can use,
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info".
% 
% As a counterpart to the access to the source code and rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty and the software's author,  the holder of the
% economic rights,  and the successive licensors have only limited
% liability.
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean that it is complicated to manipulate,  and that also
% therefore means that it is reserved for developers and experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or
% data to be ensured and,  more generally, to use and operate it in the
% same conditions as regards security.
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.

% [gradient, lipschitzConstant, convolution, convolution_T] = convolutiveGradient1D(AtA, AtY, f)
% computes the gradient of g(S) = 1 / 2 ||Y - A S W||_2^2.
% with W a convolutive transform
% Inputs:
%  - AtA: matrix A' * A.
%  - AtY: matrix A' * Y.
%  - f: filter of the convolutive transform.
% Outputs:
%  - gradient: gradient of g (handle function)
%  - lipschitzConstant: Lipschitz constant of the gradient of g.
%  - convolution: convolution handle function.
%  - convolution_T: handle function of the transpose of the convolution.

function [gradient, lipschitzConstant, convolution, convolution_T] = convolutiveGradient1D(AtA, AtY, f)
cv_F = @(F, x) real(ifft(fft(x, [], 2) .* (ones(size(x, 1), 1) * F), [], 2));
F = fft(f, [], 2);
if nargout>2
    convolution = @(x) real(ifft(fft(x, [], 2) .* (ones(size(x, 1), 1) * F), [], 2));
end
if nargout>3
    convolution_T = @(x) real(ifft(fft(x, [], 2) .* (ones(size(x, 1), 1) * conj(F)), [], 2));
end
lipschitzConstant = max(eig(AtA)) * sum(cv_F(conj(F), f));


AtYLt =cv_F(conj(F), AtY);
FFc = conj(F) .* F; clear F;

gradient = @(x)  AtA ...
    * ifft((ones(size(x, 1), 1) * FFc)...
    .* fft(x, [], 2), [], 2)...
    - AtYLt;

end



