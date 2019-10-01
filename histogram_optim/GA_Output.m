function [state,options,optchanged] = GA_Output(options,state,flag)
%GA_Output GA cusmtom output function for histogram optimization.
%   [state, options, optchanged] = GA_Output(options,state,flag)
%
%   STATE: A structure containing the following information about the state 
%   of the optimization for GA:
%             Generation: Current generation number
%              StartTime: Time when GA started, returned by TIC
%               StopFlag: Character vector containing the reason for stopping
%        LastImprovement: Generation at which the last improvement in
%                         fitness value occurred
%    LastImprovementTime: Time at which last improvement occurred
%                   Best: Vector containing the best score in each generation
%                    how: String describing the 'auglag' nonlinear
%                         constraint step
%                FunEval: Cumulative number of function evaluations
%            Expectation: Expectation for selection of individuals
%              Selection: Indices of individuals selected for elite,
%                         crossover and mutation
%             Population: Population in the current generation
%                  Score: Scores of the current population
%             NonLinIneq: Nonlinear inequality constraint violations,
%                         exists for noninteger problems with nonlinear constraints
%              NonLineEq: Nonlinear equality constraint violations,
%                         exists for noninteger problems with nonlinear constraints
%       LinearConstrType: Type of linear constraints, one of 'boundconstraints',
%                         'linearconstraints', or 'unconstrained'
%         IsMixedInteger: Boolean value, true for integer-constrained problems
%
%
%   STATE for GAMULTIOBJ:
%             Generation: Current generation number
%              StartTime: Time when GA started, returned by TIC
%               StopFlag: Character vector containing the reason for stopping
%                FunEval: Cumulative number of function evaluations
%              Selection: Indices of individuals selected for elite,
%                         crossover and mutation
%                  mIneq: Number of nonlinear inequality constraints
%                    mEq: Number of nonlinear equality constraints
%                   mAll: Number of nonlinear constraints = mIneq + mEq
%             Population: Population in the current generation
%                      C: Nonlinear inequality constraints at current point
%                    Ceq: Nonlinear equality constraints at current point
%                 isFeas: Feasibility of population, a logical vector
%           maxLinInfeas: Maximum infeasibility with respect to linear constraints
%                  Score: Scores of the current population
%                   Rank: Vector of ranks of population
%               Distance: Vector of distances of each member of the population
%                         to the nearest neighboring member
%        AverageDistance: Average of Distance
%                 Spread: Vector of spread in each generation
%       LinearConstrType: Type of linear constraints, one of 'boundconstraints',
%                         'linearconstraints', or 'unconstrained'
%         IsMixedInteger: false (gamultiobj does not support integer constraints)
%
%   FLAG: Current state in which OutputFcn is called. Possible values are:
%         init: initialization state 
%         iter: iteration state
%    interrupt: subproblem for nonlinear constraints state
%         done: final state
%
%   OPTIONS: Options object as created by OPTIMOPTIONS.
% 		
%   OPTCHANGED: Boolean indicating if the options have changed. Set to true
%   to have new the solver use new options.
%
%   Stop the solver by setting STATE.StopFlag to a nonempty value such as
%   'Y'.
%


optchanged = false;

% switch flag
%     case 'init'
%         disp('Starting the algorithm');
%     case {'iter','interrupt'}
%         disp('Iterating ...')
%     case 'done'
%         disp('Performing final task');
% end

global ga_state;
ind = length(ga_state) + 1;
ga_state{ind}.population = state.Population;
ga_state{ind}.score = state.Score;

