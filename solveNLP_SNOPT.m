function [x_opt, f_opt, information] = solveNLP_SNOPT(problem, options)

% This function is given a nonlinear program of the form
%    min f(x)  s.t. xl <=   x  <= xu
%                   bl <=  A*x <= bu
%                   cl <= c(x) <= cu
% and solves it using SNOPT.

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

% If you want to use gradient information, additionally provide an options 
% struct with the fields
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
options.solver = 'snopt';


%% check problem data for completeness and set up missing entries using default values

% problem dimensions
    % n_x     number of variables
    % n_lin   number of linear constraints
    % n_nln	  number of nonlinear constraints

[problem, n_x, n_lin, n_nln] = setupNLP_missingData(problem);


%% rewrite options into snopt form

if options.objectiveGradient == true && options.constraintsJacobian == true
    snseti('Derivative level',3)
elseif options.objectiveGradient == false && options.constraintsJacobian == true
    snseti('Derivative level',2)
elseif options.objectiveGradient == true && options.constraintsJacobian == false
    snseti('Derivative level',1)
else
    snseti('Derivative level',1)
end

% check the derivatives for errors
% snseti('Verify level',3)


%% reformulate the NLP in snopt form

% Combine objective function, nonlinear and linear constraints into one
% function:
    % F = [f(x); c(x); A*x;]
    % DF = [Df; Dc; A]

% initial value for x and box constraints on x
x_start = problem.x_start;
x_lower = problem.xl;
x_upper = problem.xu;

% lower and upper bounds for objective function and constraints
F_lower = [-inf; problem.cl; problem.bl];
F_upper = [inf; problem.cu; problem.bu];
Objective_add = 0; % additional constant term of the objective function
Objective_row = 1; % row of the objective function in F

% indices and entries of constant part of DF = [0; 0; A]
[iAfun, jAvar, A] = find([zeros(1,n_x);...
                          zeros(n_nln, n_x);...
                          problem.A]);

% indices of nonconstant part of DF = [Df; Dc; 0]
[iGfun, jGvar] = find([ones(1,n_x);...
                       ones(n_nln, n_x);...
                       zeros(n_lin, n_x)]);


% some multipliers, initialized with zero
x_multipliers = zeros(n_x,1); % initial multipliers of variables
x_state = zeros(n_x,1); % initial state of variables (basis?)
F_multipliers = zeros(1+ n_nln + n_lin, 1); % initial multipliers of constraints
F_state = zeros(1 + n_nln + n_lin, 1); %initial state of constraints (basis?)


%% call of the solver snopt
% [x,F,inform,xmul,Fmul,xstate,Fstate,output] = snopt( x, xlow, xupp, xmul, xstate,...
%                                                      Flow, Fupp, Fmul, Fstate, userfun,...
%                                                      ObjAdd, ObjRow,...
%                                                      A, iAfun, jAvar, iGfun, jGvar,...
%                                                      options)

[x,F,inform,~,~,~,~,~] = snopt (x_start, x_lower, x_upper, x_multipliers, x_state,...
                                F_lower, F_upper, F_multipliers, F_state, @objective_constraints,...
                                Objective_add, Objective_row,...
                                A, iAfun, jAvar, iGfun, jGvar );
                                                    
%% compute return values
x_opt = x(1:n_x,1);
f_opt = F(1);

information.message = inform;
information.maxVio_box = max([max(x_opt-problem.xu, 0);...
                               max(problem.xl-x_opt, 0)]);
information.maxVio_lin = max([max(problem.A*x_opt-problem.bu, 0);...
                               max(problem.bl-problem.A*x_opt, 0)]);
information.maxVio_nln = max([max(problem.nlcons(x_opt)-problem.cu, 0);...
                               max(problem.cl-problem.nlcons(x_opt), 0)]);


%% auxiliary functions

function [F, DF] = objective_constraints(x)
    % This function combines the objective function and the nonlinear
    % constraints into one function with
        % function = [f(x); c(x); A*x]
        % jacobian = [Df; Dc; A]
    % But it only evaluates the objective function and the nonlinear 
    % constraints, not the linear constraints:
        % F = [f(x); c(x); 0]
        % DF = [Df; Dc; 0]

    if nargout == 1
        % compute only the function values
        f = problem.objective(x);
        c = problem.nlcons(x);

        F = [f; c; zeros(n_lin,1)];
    else 
        % also compute derivatives, if possible
        if options.objectivegradient == true
            [f, Df] = problem.objective(x);
        else
            f = problem.objective(x);
            Df = nan(1,n_x);
        end
        
        if options.constraintsjacobian == true
            [c, Dc] = problem.nlcons(x);
        else
            c = problem.nlcons(x);
            Dc = nan(n_nln,n_x);
        end

        F = [f; c; zeros(n_lin,1)];
        DF = [Df; Dc]; 
        % sort entries of DF columnwise into one long vector, such that the
        % order coincides with iGfun, jGvar
        DF = DF(:);
    end
end

end