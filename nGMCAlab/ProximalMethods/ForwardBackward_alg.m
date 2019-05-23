% ForwardBackward_alg.m - This file is part of nGMCALab.
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
% function algorithm = ForwardBackward_alg(parameters)
% This function creates a structure triggered by the function ApplyAlgorithm.
% Inputs can be provided through a parameter structure at the creation of the
% algorithm (calling this function) or when using it (calling ApplyAlgorithm).
%
% The Forward-Backward algorithm (see Combettes and Wajs 2005) solves
% min_x f(x)+\lambda g(x) with f convex differentiable with Lipschitz gradient and
% g proper convex and lower semi-continuous function.
% It requires the following fields in the parameter structure (when calling
% ApplyAlgorithm if not set before):
% - proximal: proximal operator of g under the form @(x,lambda) prox.
% - lambda: the lambda parameter in the above problem
% - gradient: gradient of f under the form @(x) grad.
% - LipschitzConstant: Lipschitz constant of the gradient of f.
% - initialization.x: initialization of the algorithm.
%
%
%% Exampe of code to use the Forward-Backward algorithm
% % create data
% x=rand(300,1).*(mod(1:300,4)'<1);
% A=randn(200,300);
% y=A*x;
%
% %% Solve argmin_x 0.5* ||y-A x||_2^2 + lambda || x ||_1 
% 
% clear parameters
% lambda=0.5;
% %convervence properties
% parameters.MaximumIteration=500;
% parameters.RelativeDifferenceTolerance=0.0001;
% % record the function value
% parameters.recording.funcVal=@(data) 0.5*norm(al(y-A*data.x))^2+lambda*sum(abs(data.x(:)));
% % display it during the iterations
% parameters.display=@(data) plot(data.recording.funcVal);
% 
% % function information
% parameters.lambda=lambda;
% parameters.LipschitzConstant=max(eig(A'*A));
% parameters.gradient=@(x) A'*(A*x-y);
% parameters.proximal=@(x,threshold) softThresholding(x,threshold);
% parameters.initialization.x=zeros(300,1);
% 
% %launch the algorithm
% result = ApplyAlgorithm(ForwardBackward_alg(),parameters);
%
%
function algorithm = ForwardBackward_alg(parameters)

% if no parameter is given, used default ones
if nargin<1
    parameters=[];
end
algorithm.parameters=parameters;

%%
algorithm.initialize=@m_algo_initialize;
algorithm.iterate=@m_algo_iterate;
algorithm.terminate=@m_algo_terminate;
algorithm.name='Forward Backward';
end



%% ----------------------------------------------------------------------%%

function [data,parameters]=m_algo_initialize(data,parameters)
% initialization of the algorithm


%% check for required parameters
if ~isfield(parameters,'proximal')
    error('In ForwardBackward: You need to input the proximal operator in the parameter structure.\n')
end
if ~isfield(parameters,'gradient')
    error('In ForwardBackward: You need to input the gradient function in the parameter structure.\n')
end
if ~isfield(parameters,'LipschitzConstant')
    error('In ForwardBackward: You need to input LipschitzConstant for the gradient in the parameter structure (2*lambda_max).\n')
end
if ~isfield(data,'x')
    error('In ForwardBackward: You need to input initialization.x in the parameter structure.\n')
end

%% initialize automatically missing general parameters
parameters = m_automaticSettings('lambda',1,parameters);
parameters = m_automaticSettings('MaximumIteration',10000,parameters);
parameters = m_automaticSettings('RelativeDifferenceTolerance',10^-6,parameters);


%% initialization
data.previous_x=data.x;
data.t=1;

end


%% ----------------------------------------------------------------------%%

function [data,relativeDifference]=m_algo_iterate(data,parameters)
% iterate and update the data (and optionnally the function value)

    data.t_next=(1+sqrt(1+4*data.t^2))/2;
    w=(data.t-1)/data.t_next;
    data.t=data.t_next;
    y=(1+w)*data.x-w*data.previous_x;

data.previous_x=data.x;
data.x=parameters.proximal(...
    y-parameters.gradient(y)/parameters.LipschitzConstant,...
    parameters.lambda/parameters.LipschitzConstant);


%data.t=(1+sqrt(1+4*data.t^2))/2;

%%
relativeDifference=norm(al(data.x-data.previous_x))/norm(al(data.previous_x));

end






%% ----------------------------------------------------------------------%%

function result=m_algo_terminate(data)
% output the result of the algorithm

result.x=data.x;
end




%% ------------------------- INTERNAL FUNCTIONS -------------------------%%
%% ----------------------------------------------------------------------%%

function parameters = m_automaticSettings(name,value,parameters)
if ~isfield(parameters,name)
    parameters.(name)=value;
    if isfield(parameters,'display')
        fprintf(1,'In %s, parameter "%s" automatically set to %f.\n',parameters.name,name,parameters.(name));
    end
end
end



