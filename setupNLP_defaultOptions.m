function options = setupNLP_defaultOptions(options)

% This function is given an empty argument or an options struct with some
% or all of the following fields:
    % options.objectiveGradient
    % options.constraintsJacobian
    % options.solver
    
% It sets up missing or empty fields using the following default values:
    % options.objectiveGradient = false
    % options.constraintsJacobian = false
    % options.solver = 'fmincon'

%%

if isempty(options)
    options.objectiveGradient = [];
    options.constraintsJacobian = [];
    options.solver = [];
end

if ~isfield(options, 'objectiveGradient') || isempty(options.objectiveGradient)
    % default value is no gradient information
    options.objectiveGradient = false;
end

if ~isfield(options, 'constraintsJacobian') || isempty(options.constraintsJacobian)
    % default value is no gradient information
    options.constraintsJacobian = false;
end

if ~isfield(options, 'solver') || isempty(options.solver)
    % default solver is fmincon
    options.solver = 'fmincon';
end