% analysis_nGMCA_RW_alg.m - This file is part of nGMCALab
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
% Linear operator
% Implements operator-like operations for a given function and its
% transpose
%
% M = LinearOperator(direct_function, transpose_function)
% Inputs:
% - direct_function: function handle of the direct operation
% - transpose_function: function handle of the transpose operation
%
% Output
% - M: class instance such that M * x = direct_function(x)
%      and M' * x = transpose_function(x)
% It is possible to use multiplication with a vector, a scalar, division by
% a scalar, transpose, product of operators, sum of operators just as
% though M were a matrix.
%
% With A a matrix, B = LinearOperator(@(x) A * x, @(x) A' * x);
% behaves similarly as A, except for the addition with a scalar, for which
% the scalar is considered as a identity times the scalar matrix.
classdef LinearOperator
    
    properties
        d_f;
        t_f;
    end
    
    methods
        function ope = LinearOperator(direct_function, transpose_function)
            % instanciate an operator
            ope.d_f = direct_function;
            ope.t_f = transpose_function;
        end
        
        function res = plus(a, b)
            % additions
            if isa(a, 'LinearOperator')
                % of an operator
                if isa(b, 'LinearOperator')
                    % with an operator
                    res = LinearOperator(@(x) a.d_f(x) + b.d_f(x),...
                        @(x) a.t_f(x) + b.t_f(x));
                else if isscalar(b)
                        % with a scalar (considered as scalar times the
                        % identity)
                        res = LinearOperator(@(x) a.d_f(x) + b * x,...
                            @(x) a.t_f(x) + b * x);
                    else
                        error('This plus operation with a LinearOperator is not implemented.\n')
                    end
                end
            else
                if isscalar(a)
                    % of a scalar with an operator(considered as scalar
                    % times the identity)
                    res = LinearOperator(@(x) b.d_f(x) + a * x,...
                        @(x) b.t_f(x) + a * x);
                else
                    error('This plus operation with a LinearOperator is not implemented.\n')
                end
            end
            
        end
        
        
        
        function res = minus(a, b)
            % substration
            if isa(a, 'LinearOperator')
                % of an operator
                if isa(b, 'LinearOperator')
                    % with an operator
                    res = LinearOperator(@(x) a.d_f(x) - b.d_f(x),...
                        @(x) a.t_f(x) - b.t_f(x));
                else if isscalar(b)
                        % with a scalar
                        res = LinearOperator(@(x) a.d_f(x) - b * x,...
                            @(x) a.t_f(x) - b * x);
                    else
                        error('This minus operation with a LinearOperator is not implemented.\n')
                    end
                end
            else
                if isscalar(a)
                    % of a scalar with an operator
                    res = LinearOperator(@(x) -b.d_f(x) + a * x,...
                        @(x) -b.t_f(x) + a * x);
                else
                    error('This minus operation with a LinearOperator is not implemented.\n')
                end
            end
            
        end
        
        
        
        function res = mtimes(a, b)
            % product
            if isa(a, 'LinearOperator')
                % of an operator
                if isa(b, 'LinearOperator')
                    % with an operator (concatenation of operator)
                    res = LinearOperator(@(x) a.d_f(b.d_f(x)),...
                        @(x) b.t_f(a.t_f(x)));
                else if isscalar(b)
                        % with a scalar
                        res=LinearOperator(@(x) b * a.d_f(x),...
                            @(x)  b * a.t_f(x));
                    else
                        % with an object such as vectors or matrices
                        % then apply the direct function on it
                        res = a.d_f(b);
                    end
                end
            else
                if isscalar(a)
                    % of a scalar with an operator
                    res = LinearOperator(@(x) a * b.d_f(x),...
                        @(x)  a * b.t_f(x));
                else
                    % of an object x with an operator M
                    % compute compute x * M = (M' * x')'
                    res = b.t_f(a')';
                end
                
            end
            
        end
        
        function res = mrdivide(a, b)
            % division
            if isa(a, 'LinearOperator') && isscalar(b)
                % of an operator by a scalar
                res = LinearOperator(@(x) a.d_f(x) / b,...
                    @(x)  a.t_f(x) / b);
            else
                error('This division operation with a LinearOperator is not implemented.\n')
            end
        end
        
        
        function res = uminus(a)
            % unary minus operator
            res = LinearOperator(-a.d_f, -a.t_f);
        end
        
        function res = ctranspose(a)
            % transpose operator
            % transpose operator swaps direct and transpose functions
            % multiplying the output operator with an object x will then 
            % compute transpose_function(x)
            res=LinearOperator(a.t_f, a.d_f);
        end
        
    end
end


