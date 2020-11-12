function [x_opt, f_opt, information] = solveNLP_FMINCON(problem, options)

% This function is given a nonlinear program of the form
%    min f(x)  s.t. xl <=   x  <= xu
%                   bl <=  A*x <= bu
%                   cl <= c(x) <= cu
% and solves it using FMINCON.

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

% If you want  to use gradient information, additionally provide an options struct 
% with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    
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
options.solver = 'fmincon';


%% check problem data for completeness and set up missing entries using default values

problem = setupNLP_missingData(problem);


%% rewrite options into fmincon form

newproblem.options = [];

% specify, if gradient information is provided for the objective
if options.objectiveGradient == true
    newproblem.options.SpecifyObjectiveGradient = true;
else 
    newproblem.options.SpecifyObjectiveGradient = false;
end

% specify, if gradient information is provided for the constraints
if options.constraintsJacobian == true
    newproblem.options.SpecifyConstraintGradient = true;
else
    newproblem.options.SpecifyConstraintGradient = false;
end

% have Matlab check the gradient for correctness
if (options.objectiveGradient == true) || (options.constraintsJacobian == true)
    newproblem.options.CheckGradients = true;
end

% supress output from fmincon
newproblem.options.Display = 'off';

% choose the solution algorithm used by fmincon
newproblem.options.Algorithm = 'sqp';

% increase iteration limits for fmincon
newproblem.options.MaxFunctionEvaluations = 10^7;
newproblem.options.MaxIterations = 10^4;



%% reformulate the NLP in fmincon form

% objective function
newproblem.objective = problem.objective;

% initial value
newproblem.x0 = problem.x_start;

% box constraints
newproblem.lb = problem.xl;
newproblem.ub = problem.xu;

% linear inequality and equality constraints, need to be separated for
% FMINCON into
    % Aineq * x <= bineq
    % Aeq   * x  = beq
eq_lin = (problem.bl == problem.bu);
ineq_lin_upper = logical((problem.bu <  inf).*(problem.bu ~= problem.bl));
ineq_lin_lower = logical((problem.bl > -inf).*(problem.bl ~= problem.bu));
newproblem.Aineq = [  problem.A(ineq_lin_upper,:); ...
                    - problem.A(ineq_lin_lower,:)];
newproblem.bineq = [  problem.bu(ineq_lin_upper); ...
                    - problem.bl(ineq_lin_lower)];
newproblem.Aeq = problem.A(eq_lin,:);
newproblem.beq = problem.bl(eq_lin);

% nonlinear inequality and equality constraints
newproblem.nonlcon = @nonlinearConstraints;

newproblem.solver = 'fmincon';


%% solve the NLP with fmincon

[x_opt, f_opt, ~, output] = fmincon(newproblem);


%% compute return values

information.message = output.message;
information.maxVio_box = max([max(x_opt-problem.xu, 0);...
                               max(problem.xl-x_opt, 0)]);
information.maxVio_lin = max([max(problem.A*x_opt-problem.bu, 0);...
                               max(problem.bl-problem.A*x_opt, 0)]);
information.maxVio_nln = max([max(problem.nlcons(x_opt)-problem.cu, 0);...
                               max(problem.cl-problem.nlcons(x_opt), 0)]);

%% auxiliary functions

function [cineq, ceq, Dcineq, Dceq] = nonlinearConstraints(x)
    % need to be separated for FMINCON into
        % cineq(x) <= 0
        % ceq(x)    = 0
    
    if nargout == 2
        % only functions evaluated
        c = problem.nlcons(x);
        
        eq_nln = (problem.cl == problem.cu);
        ineq_nln_upper = logical((problem.cu <  inf).*(problem.cu ~= problem.cl));
        ineq_nln_lower = logical((problem.cl > -inf).*(problem.cl ~= problem.cu));
        
        ceq = c(eq_nln) - problem.cl(eq_nln);
        cineq = [  c(ineq_nln_upper) - problem.cu(ineq_nln_upper);...
                 - c(ineq_nln_lower) + problem.cl(ineq_nln_lower)];

    else 
        % gradients also evaluated
        [c, Dc] = problem.nlcons(x);
        
        eq_nln = (problem.cl == problem.cu);
        ineq_nln_upper = logical((problem.cu <  inf).*(problem.cu ~= problem.cl));
        ineq_nln_lower = logical((problem.cl > -inf).*(problem.cl ~= problem.cu));
        
        ceq = c(eq_nln) - problem.cl(eq_nln);
        cineq = [  c(ineq_nln_upper) - problem.cu(ineq_nln_upper);...
                 - c(ineq_nln_lower) + problem.cl(ineq_nln_lower)];
        
        Dceq = Dc(eq_nln, :);
        Dcineq = [  Dc(ineq_nln_upper,:);...
                  - Dc(ineq_nln_lower,:)];
             
        % rotate gradients of constraints columnwise for fmincon   
        Dceq = Dceq';
        Dcineq = Dcineq';
    end
end
    

end