% LoopTimer.m - This file is part of nGMCALab
% This software aims at performing non-negative matrix factorization.
% Copyright 2013 CEA
% Contributor : Jérémy Rapin (jeremy.rapin@cea.fr)
% Created on 20/08/13
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
% LoopTimer class
% Used in performBenchmark.m in order to evaluate the remaining running time.
% BEWARE: buggy but good enough most of the time.
%
classdef LoopTimer < handle
    properties
        times
        iteration
        n_iterations
        clockID
    end
    
    methods
        
        function timer=LoopTimer(number_of_iterations)
            timer.n_iterations=number_of_iterations;
            timer.iteration=0;
            timer.times=NaN*ones(number_of_iterations,1);
            timer.clockID=0;
        end
        
        
        
        function notifyNewIteration(timer)
            
            if timer.iteration==0
                timer.clockID=tic();
            end
            
            timer.iteration=timer.iteration+1;
            timer.times(timer.iteration)=toc(timer.clockID);
        end
        
        function secs=secondsRemaining(timer)
            curr=timer.iteration;
            last=timer.n_iterations;
            firs=max(1,curr-(last+1-curr));%first time used for approximation
            %at the very end, no point using the times from the very beginning
            
            if curr>1
                speed=(timer.times(curr)-timer.times(firs))/(curr-firs);
                secs=(last+1-curr)*speed-(toc(timer.clockID)-timer.times(curr));
            else
                secs=NaN;
            end
            
        end
        
        
        function str=stringRemainingTime(timer)
            
            secs=timer.secondsRemaining();
            
            if isnan(secs)
                str='unknown';
            else
                seconds=mod(floor(secs),60);
                str=[num2str(seconds) 's'];
                minutes=mod(floor(secs/60),60);
                if minutes>0
                    str=[num2str(minutes) 'min' str];
                    hours=mod(floor(secs/3600),24);
                    if hours>0
                        str=[num2str(hours) 'h' str];
                        days=floor(secs/86400);
                        if days>0
                            str=[num2str(days) 'd' str];
                        end
                    end
                end
            end
        end
    end
end

