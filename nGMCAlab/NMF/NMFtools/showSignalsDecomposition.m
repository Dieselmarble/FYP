% showSignalsDecomposition.m - This file is part of nGMCALab.
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13, last modified on 16/7/2014
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software. You can use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and rights to copy, 
% modify and redistribute granted by the license, users are provided only
% with a limited warranty and the software's author, the holder of the
% economic rights, and the successive licensors have only limited
% liability. 
%
% In this respect, the user's attention is drawn to the risks associated
% with loading, using, modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software, 
% that may mean that it is complicated to manipulate, and that also
% therefore means that it is reserved for developers and experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and, more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.
%
%
%
% showImageDecomposition(result, reference)
% Shows the decomposition of the first 4 signals of results.S in interference, 
% noise and artifacts as defined in Vincent et al., 2006
% can be used as display function in ApplyAlgorithm by setting:
% parameters.display = @(data) showImageDecomposition(data, reference)
%
% Inputs:
% - result: the structure of the results with fields 'A' and 'S'
% - reference: structure containing reference 'A' and 'S'.


function showSignalsDecomposition(data, reference)
data = NMFnormalization(data, 'S');
[criteria, data] = NMFevaluation(data, reference, 0);
 d = data.decomposition;
 clf
 
 subplot(2, 2, 1)
 m_decompPlot(d, 1, 0)
 subplot(2, 2, 2)
 m_decompPlot(d, 2, 0)
 subplot(2, 2, 3)
 m_decompPlot(d, 3, 0)
 subplot(2, 2, 4)
 m_decompPlot(d, 4, 1)


end


function m_decompPlot(d, k, leg)
 hold on
 
 M = max([d.target(k, :)'; d.artifacts(k, :)'; d.noise(k, :)'; d.interferences(k, :)';...
 d.artifacts(k, :)' + d.noise(k, :)' + d.interferences(k, :)' + d.target(k, :)']);
 
 plot(d.target(k, :)' / M, 'color', 'blue', 'LineWidth', 2);
 plot(d.artifacts(k, :)' / M, 'color', 'cyan', 'LineWidth', 2);
 plot(d.noise(k, :)'/M, 'color', 'green', 'LineWidth', 2);
 plot(d.interferences(k, :)'/M, 'color', 'red', 'LineWidth', 2);
 plot((d.artifacts(k, :)'+d.noise(k, :)'+d.interferences(k, :)'+d.target(k, :)')/M, '--', 'color', 'black', 'LineWidth', 2);
 v = axis();
 v(2) = size(d.artifacts, 2);
 v(4) = 1.05;
 axis(v)
 
 if leg
 LEG = legend('target', 'artifacts', 'noise', 'interferences', 'reconstruction');
 set(LEG, 'FontSize', 13, 'FontName', 'Times New Roman');
 end
 hold off
 
 title(['Source ' num2str(k)], 'FontSize', 13, 'FontName', 'Times New Roman');
 
 set(gca, ...
 'Box', 'off' , ...
 'TickDir', 'out' , ...
 'XMinorTick', 'on' , ...
 'YMinorTick', 'off' , ...
 'YGrid', 'off' , ...
 'YTick', 0 , ...
 'XColor', [.15 .15 .15], ...
 'YColor', [.15 .15 .15], ...
 'LineWidth', 1, ...
 'FontName', 'Times New Roman', ...
 'FontSize', 13);
end


