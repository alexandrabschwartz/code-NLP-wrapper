# code-NLP-wrapper
Wrapper to provide a unified call for different NLP solvers in MATLAB. The included functions are described below. 

I did test this code, but it may still contain bugs. If you find errors or have suggestions for improvements, please let me know.

## solveNLP
This function takes a nonlinear optimization problem as well as optional options as input. Within the options, you can specify which NLP solver you want to use. At them moment, FMINCON, FMINUNC and SNOPT are possible. The function then passes the problem on to a call-function for the specified NLP solver.

## solveNLP_FMINCON
This function takes a nonlinear optimization problem as well as optional options as input. It rewrites the optimization problem and the options into FMINCON form and then calls FMINCON to solve the NLP.

## solveNLP_FMINUNC
This function takes a nonlinear optimization problem as well as optional options as input. It rewrites the optimization problem and the options into FMINUNC form and then calls FMINUNC to solve the NLP. FMINUNC is a solver for unconstrained problem, so you shoudl not use thos sovler if you optimization problem has constraints.

## solveNLP_SNOPT
This function takes a nonlinear optimization problem as well as optional options as input. It rewrites the optimization problem and the options into SNOPT form and then calls SNOPT to solve the NLP. SNOPT is a paid-for NLP solver, but limited trial versions can be obtained for free. This function has not been tested extensively.

## setupNLP_missingData
This function takes a nonlinear optimization problem as input. It checks thee problem data for completeness and inserts missing data -- if possible -- using default values. E.g. if you did not specify box constraints on the variable it inserts -inf/inf as lower/upper bounds.

## setupNLP_defaultOptions
This function takes an options struct as input and sets up missing options using default values. E.g. if you did not specify the NLP solver, it chooses FMINCON.

## testNLP
This script contains two toy nonlinear problems to illustrate how the functions are called.
