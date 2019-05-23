% External file
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
%  Get NMR spectra from the SDBS database peak info (http://riodb01.ibase.aist.go.jp/sdbs/)
%
%  function [y,x] = Get_NMR_Spectrum(peak,ppm,prange,nsamples,fwhm)
%  
%  peak : peak amplitude
%  ppm : value of the ppm (position of the peaks)
%  prange = [pmin,pmax] : desired range of ppm
%  nsamples : desired number of samples
%  fwhm : full width at half maximum of the Laplacian kernel (profile of the peaks)
%
%
%   Example :  H1-NMR for the mannitol :
%
%	ppm=[4.432  4.418  4.362  4.347  4.333  4.163  4.145  3.637  3.629  3.623  3.614  3.610  3.602  3.596  3.587  3.562  3.543  3.524  3.485  3.476  3.470  3.462  3.456  3.449  3.441  3.435  3.427  3.408  3.394  3.381  3.366  3.352];
%peak = [  965 1000  382  917  409  772  819  164 191  172  225  236  224  215  217  236  470  357  104  112  212  212 198  217  128  110  101  263  343  275  299 148];
%
%
%   [y,x] = Get_NMR_Spectrum(peak,ppm,[1,10],10000,5);
%
%
%
%


function [y,x] = Get_NMR_Spectrum(peak,ppm,prange,nsamples,fwhm)

%----

x = (max(prange)-min(prange))*(0:nsamples-1)/(nsamples-1) + min(prange);
y = 0.*x;

for npeak = 1:numel(peak)
	
	ii = find(abs(x-ppm(npeak)) == min(abs(x-ppm(npeak))));
	ii = ii(1);
	y(ii) = y(ii) + peak(npeak);
	
end

%----


np = floor(nsamples/2);
kern = [np:-1:1,0,1:np];kern = kern(1:nsamples);
lambda = 0.5*fwhm/log(2);
kern = exp(-kern/lambda);

y = conv(y,kern);%,'same');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if length(y)>length(kern)
   L= floor(length(kern)/2.0);
   y=y(L+1:end-L+1);
end

end


