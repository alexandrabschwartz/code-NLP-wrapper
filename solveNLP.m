function [x_opt, f_opt, information] = solveNLP(problem,options)

% This function is given a nonlinear program of the form
%    min f(x)  s.t. xl <=   x  <= xu
%                   bl <=  A*x <= bu
%                   cl <= c(x) <= cu
% and solves it using FMINCON or another NLP solver.

% The nonlinear program should be provided as a struct with the following
% fields: 
    % problem.objective = @objective;
    % problem.xl = xl;
    % problem.xu = xu; 
    % problem.A = A;
    % problem.bl = bl;
    % problem.bu = bu;
    % problem.nlcons = @nlcons;
    % problem.cl = cl;
    % problem.cu = cu;
    % problem.x_start = x_start;
    % problem.dimension = n_x ;
% For the objective function and the nonlinear constraint the respective 
% functions can either only return the function value or additionally the 
% gradients (oriented row-wise). The default assumption is that no gradients
% are provided.

% If you want to specify the NLP solver or use gradient information, additionally 
% provide an options struct with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.solver = 'fmincon' or 'fminunc' or 'snopt'
    
% The function returns
    % x_opt                    computed solution
    % f_opt                    objective function value in x_opt
    % information.message      exit message of the solver
    % information.maxVio_box   maximum violation of box constraints
    % information.maxVio_lin   maximum violation of linear constraints
    % information.maxVio_nln   maximum violation of nonlinear constraints


%% set up options using default values

if nargin == 1
    options = [];
end
options = setupNLP_defaultOptions(options);



%% call the specified NLP solver
switch options.solver
    case 'fmincon'
        [x_opt, f_opt, information] = solveNLP_FMINCON(problem, options);
    case 'snopt'
        [x_opt, f_opt, information] = solveNLP_SNOPT(problem, options);
    case 'fminunc'
        [x_opt, f_opt, information] = solveNLP_FMINUNC(problem, options);
    otherwise
        disp('NLP solver unknown, using fmincon instead')
        [x_opt, f_opt, information] = solveNLP_FMINCON(problem, options);
end
    