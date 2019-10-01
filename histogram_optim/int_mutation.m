%---------------------------------------------------
% Mutation function to generate childrens satisfying the range and integer
% constraints on decision variables.
function mutationChildren = int_mutation(parents,options,GenomeLength, ...
    FitnessFcn,state,thisScore,thisPopulation)
shrink = .01; 
scale = 1;
scale = scale - shrink * scale * state.Generation/options.Generations;
range = options.PopInitRange;
lower = range(1,:);
upper = range(2,:);
scale = scale * (upper - lower);
mutationPop =  length(parents);
% The use of ROUND function will make sure that childrens are integers.
mutationChildren =  repmat(lower,mutationPop,1) +  ...
    round(repmat(scale,mutationPop,1) .* rand(mutationPop,GenomeLength));
% End of mutation function
%---------------------------------------------------