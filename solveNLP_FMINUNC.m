function [x_opt, f_opt, information] = solveNLP_FMINUNC(problem, options)

% This function is given an unconstrained nonlinear program of the form
%    min f(x)  s.t. x in R^n
% and solves it using FMINUNC.

% The unconstrained nonlinear program should be provided as a struct with
% the following fields: 
    % problem.objective = @objective;
    % problem.x_start = x_start;
    % problem.dimension = n_x ;
% The objective function can either only return the function value or
% additionally the gradient (oriented row-wise).
% The default assumption is that no gradient is provided.

% If you want  to use gradient information, additionally provide an options
% struct with the field
    % options.objectiveGradient = true or false
    
% The function returns
    % x_opt                    computed solution
    % f_opt                    objective function value in x_opt
    % information.message      exit message of the solver


%% set up options using default values

if nargin == 1
    options = [];
end
options = setupNLP_defaultOptions(options);
options.solver = 'fminunc';


%% check problem data for completeness and set up missing entries using default values

problem = setupNLP_missingData(problem);


%% rewrite options into fminunc form

newproblem.options = [];

% specify, if gradient information is provided for the objective
if options.objectiveGradient == true
    newproblem.options.SpecifyObjectiveGradient = true;
else 
    newproblem.options.SpecifyObjectiveGradient = false;
end

% have Matlab check the gradient for correctness
if options.objectiveGradient == true
    newproblem.options.CheckGradients = true;
end

% supress output from fminunc
newproblem.options.Display = 'off';

% choose the solution algorithm used by fminunc
newproblem.options.Algorithm = 'quasi-newton';


%% reformulate the NLP in fminunc form

newproblem.objective = problem.objective;

newproblem.x0 = problem.x_start;

newproblem.solver = 'fminunc';


%% solve the NLP with fminunc

[x_opt, f_opt, ~, output] = fminunc(newproblem);


%% compute return values

information.message = output.message;

end