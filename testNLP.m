%%  MPVC test problem for NLP solvere
%   min 4x_1 + 2x_2 ST x >= 0,
%                       (5*sqrt(2) - x_1 - x_2)*x_1 <= 0, 
%                       (5 - x_1 - x_2)*x_2 <= 0

% [0; 0] is the global minimum and an isolated feasible point
% [0; 5] is a local minimum
% [0; 5*sqrt(2)] is not a local solution

problem.objective = @objective_MPVC;
problem.xl = [0; 0];
problem.xu = [inf; inf]; 
problem.A =[];
problem.bl = [];
problem.bu = [];
problem.nlcons = @nlcons_MPVC;
problem.cl = [-inf; -inf];
problem.cu = [0; 0];
problem.x_start = [3; 3];
problem.dimension = 2 ;

options.objectiveGradient = true;
options.constraintsJacobian = true;
options.solver = 'fmincon';

% call of the NLP solver

[x_opt, f_opt, information] = solveNLP(problem,options)




%% MPCC teestproblem for NLP solver
%   min (x_1-1)^2 + (x_2-1)^2 ST x >= 0,
%                                x_1 * x_2 <= 0

% global minima are [1; 0] and [0; 1]

problem.objective = @objective_MPCC;
problem.xl = [0; 0];
problem.xu = [inf; inf]; 
problem.A =[];
problem.bl = [];
problem.bu = [];
problem.nlcons = @nlcons_MPCC;
problem.cl = -inf;
problem.cu = 0;
problem.x_start = [1; 0];
problem.dimension = 2 ;

options.objectiveGradient = true;
options.constraintsJacobian = true;
options.solver = 'fmincon';


% call of the NLP solver

[x_opt, f_opt, information] = solveNLP(problem,options)



%% objective functions and nonlinear constraints for test problems

function [f_x, Df_x] = objective_MPVC(x)
    f_x = [4;2]'*x;
    
    if nargout > 1
        % gradients should be orientteed row-wise
        Df_x = [4 2];
    end
end

function [c_x, Dc_x] = nlcons_MPVC(x)
    c_x = [(5*sqrt(2)-x(1)-x(2))*x(1);...
           (5-x(1)-x(2))*x(2)];
       
    if nargout > 1
        % gradients should be orientteed row-wise
        Dc_x = [-2*x(1)+5*sqrt(2)-x(2)  -x(1);...
                 -x(2)                  -2*x(2) + 5 - x(1)];
    end
end

function [f_x, Df_x] = objective_MPCC(x)
    f_x = (x(1)-1)^2 + (x(2)-1)^2;
    
    if nargout > 1
        % gradients should be orientteed row-wise
        Df_x = 2*[x(1)-1  x(2)-1];
    end
end

function [c_x, Dc_x] = nlcons_MPCC(x)
    c_x = x(1) * x(2);
    
    if nargout > 1
        % gradients should be orientteed row-wise
        Dc_x = [x(2)  x(1)];
    end
end

        