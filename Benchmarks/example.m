clearvars
close all
clc

% This file has examples of solving the two provided benchmark problems
% using MATLAB's "surrogateopt" and "fmincon".  It can be used as a model 
% for how to solve the problems with other optimizers.

% Get the automotive problem

[automotive,frel,grel,hrel,xl,xu,xopt,x0] = automotive_benchmark();

% The output "automotive" is a function that is called as follows:
%
%     [f,g,h] = automotive(x)
%
% where
%
%     x = row vector of size 1 x d where d = length(xl).
%         The vector x should obey the bounds xl <= x <= xu.
%     f = scalar objective to be minimized.
%     g = vector of inequality constraints of size 1 x ng,
%         where ng = size(grel,1).  We require g <= 0.
%     h = vector of equality constraints of size 1 x nh,
%         where nh = size(hrel,1).  We require h = 0.
%
%         NOTE: for the automotive problem, hrel=[]; there are no equalities.
%
% The other outputs of 'automotive_benchmark' are:
%
%     frel = indicates which variables have an impact on the
%            objective.  If frel(k)=true, then variable k affects f;
%            if frel(k)=false, then variable k has no impact on f.
%
%     grel = indicates which variables have an impact on the inequalities.
%            If grel(i,k)=true, then variable k affects g(i);
%            if grel(i,k)=false, then variable k has no impact on g(i).
%
%     hrel = indicates which variables have an impact on the equalities.
%            If hrel(i,k)=true, then variable k affects h(i);
%            if hrel(i,k)=false, then variable k has no impact on h(i).
%
%            NOTE:  The information on which variables affect the different
%            outputs can be used when fitting surrogate models; 
%            in particular, the non-relevant variables can be ignored when 
%            fitting the surrogate.
%
%     xl   = vector of size 1 x d with lower bounds on x;
%
%     xu   = vector of size 1 x d with upper bounds on x;
%
%     xopt = presumed global minimum

% Setup the problem for solving with Matlab's surrogateopt.

objconstr = @(x) objconstr_for_problem(x,automotive);

options = optimoptions('surrogateopt',Display='iter',...
                       MaxFunctionEvaluations=4000,...
                       InitialPoints=x0);

% Set a random number seed for repeatability.
% Note that surrogateopt's performance is sensitive to the seed.

rng(43);

% Call surrogateopt

xmin = surrogateopt(objconstr,xl,xu,options);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the business jet problem.
% Use the IDF version that HAS equality constraints.
% Use sigmoid scaling for the constraints.

UseEq = true;
Scale = true;

[bjet,frel,grel,hrel,xl,xu,xopt] = businessjet_benchmark(UseEq,Scale);

% The problem function 'bjet' is called as follows:
%
%     [f,g,h] = bjet(x)
%
% where the meanings of {x,f,g,h,frel,grel,hrel,xl,xu,xopt} are the same as
% described for the automotive problem.

% Set up the problem for solving with fmincon

fun     = @(x) fun_for_problem(x,bjet);
nonlcon = @(x) nonlcon_for_problem(x,bjet);

options = optimoptions('fmincon',Algorithm='sqp',...
                       FiniteDifferenceType='central',Display='iter',...
                       MaxFunctionEvaluations=83000,MaxIterations=1000);

% Load a saved starting point  

load('x0.mat','x0');

% From this point, fmincon successfully converges to the global minimum.  
% But if you experiment with different random starting points, you will 
% find that the optimization very often fails to find a feasible point or
% converges to a suboptimal local minimum.

% Solve the problem with fmincon

xmin = fmincon(fun,x0,[],[],[],[],xl,xu,nonlcon,options);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the business jet problem again.
% This time use the MDF version that has NO equality constraints.
% Continue to use sigmoid scaling for the constraints.
% Note that the MDF version sometimes returns NaN for f and g when it
% fails to solve for the values of the coupling variables via functional
% iteration.  The optimization algorithm needs to be able to handle such 
% failed calculations.  In the literature, this is sometimes called a
% "hidden constraint," as there is implicitly a constraint which says that
% the calculations must be successful.

UseEq = false;
Scale = true;

[bjet,frel,grel,hrel,xl,xu,xopt] = businessjet_benchmark(UseEq,Scale);

% Solve the problem with Matlab's surrogateopt which appears to be able to 
% handle the hidden constraint. 

objconstr = @(x) objconstr_for_problem(x,bjet);

options = optimoptions('surrogateopt',Display='iter',MaxFunctionEvaluations=4000);

% Set a random number seed for repeatability.
% Note that surrogateopt's performance is sensitive to the seed.
% With the seed used here (43), it converges to a solution with f=33.47,
% which is fairly close to the global minimum of f=32.635.  But I have  
% found that other seeds give substantially worse solutions (e.g., f > 40).

rng(43);

% Call surrogateopt

xmin = surrogateopt(objconstr,xl,xu,options);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = objconstr_for_problem(x,problem)
    [f,g]  = problem(x);
    s.Fval = f;
    s.Ineq = g;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = fun_for_problem(x,problem)
    [f,~,~] = problem(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,Ceq] = nonlcon_for_problem(x,problem)
    [~,C,Ceq] = problem(x);
end
