% SDRPlot.m - This file is part of nGMCALab.
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
% SDRPlot(data)
% Plots the evolution of the criteria with the SDR, SIR, SNR, SAR during the
% iterations, as defined in Vincent et al., 2006.
% It requires recording of the criteria during the algorithm using
% ApplyAlgorithm, with parameter:
% parameters.recording.crit = @(data) NMFevaluation(data, reference);
% can be used as display function in ApplyAlgorithm by setting:
% parameters.display = @(data) SDRPlot(data);



function SDRPlot(data)
fontSize=13;

clf
subplot(2, 1, 1)
hold on
plot(data.recording.SDR_S, 'color', 'blue', 'LineWidth', 2);
plot(data.recording.SIR_S, 'color', 'red', 'LineWidth', 2);
plot(data.recording.SAR_S, 'color', 'black', 'LineWidth', 2);
plot(data.recording.SNR_S, 'color', 'green', 'LineWidth', 2);
LEG=legend('SDR_S' ,'SIR_S', 'SAR_S', 'SNR_S', 0);
set(LEG,'FontSize', fontSize, 'FontName', 'Times New Roman');

hold off
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...%'TickLength'  , [1 1]) , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.15 .15 .15], ...
    'YColor'      , [.15 .15 .15], ...
    'LineWidth'   , 1,...
    'FontSize',fontSize,...
    'FontName','Times New Roman');

xlabel('iteration', 'FontSize', fontSize, 'FontName', 'Times New Roman')
ylabel('dB', 'FontSize', fontSize, 'FontName', 'Times New Roman')

subplot(2, 1, 2)
hold on
plot(data.recording.SDR_A, 'color', 'blue', 'LineWidth', 2);
plot(data.recording.SIR_A, 'color', 'red', 'LineWidth', 2);
plot(data.recording.SAR_A, 'color', 'black', 'LineWidth', 2);
plot(data.recording.SNR_A, 'color', 'green', 'LineWidth', 2);
LEG=legend('SDR_A', 'SIR_A', 'SAR_A', 'SNR_A', 0);
set(LEG, 'FontSize', fontSize, 'FontName', 'Times New Roman');

hold off
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...%'TickLength'  , [1 1]) , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.15 .15 .15], ...
    'YColor'      , [.15 .15 .15], ...
    'LineWidth'   , 1,...
    'FontSize',fontSize,...
    'FontName','Times New Roman');

xlabel('iteration', 'FontSize', fontSize, 'FontName', 'Times New Roman')
ylabel('dB', 'FontSize', fontSize, 'FontName', 'Times New Roman')

end



