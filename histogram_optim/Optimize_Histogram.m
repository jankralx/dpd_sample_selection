function [histogram, optim_nmse] = Optimize_Histogram( pars, arch )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% stpars = pars;
% starch = arch;
% 
% st_ind = 1;

% % run architecture to find starting coefficients for histogram optimization
% starch.us.num_samples = 2 * starch.us.num_samples;
% if starch.us.num_samples < 60
%     starch.us.num_samples = 60;
% end
% 
% starch.us.hist.method = 'even';
% starch = Run_DPD_Archs(stpars, starch);
% 
% % select appropriate starting DPD coefficients
% 
% nmserr = [starch.res(:).nmse].';
% num_avg_samples = 20;
% 
% avg_nmse = mean(nmserr(end-num_avg_samples:end));
% 
% % select index of NMSE where avg_nmse is 3 dB
% st_ind = find(nmserr > avg_nmse + 3, 1, 'last');
% 
% if isempty(st_ind)
%     st_ind = 15;
% end

% arch.coefs(1,:) = starch.coefs(st_ind+1,:);     % save coefs as starting for optimizatiion


%% GA optimization

% prepare ga_state which will containt all results
global ga_state;
ga_state = cell(0);

fitness_fcn = @(x) HistogramFitness(x, pars, arch);
num_vars = 10;

lower_bound = zeros(num_vars,1);   
upper_bound = arch.us.num_samples * ones(num_vars,1);
   
bounds = [lower_bound.'; upper_bound.']; % If unbounded then Bound = []

options = gaoptimset('CreationFcn',@int_pop, 'MutationFcn',@int_mutation, ...
    'OutputFcn',@GA_Output, 'PopInitRange',bounds, 'Display','iter', ...
    'StallGenL',40, 'Generations', pars.histopt.generations, 'PopulationSize', pars.histopt.popsize, 'UseParallel', true);

[histogram,fval] = ga(fitness_fcn,num_vars,options);
fprintf('GA final histogram [');
fprintf('%d, ', histogram);
fprintf('] with NMSE = %f\n', fval);

% get global optimal histogram that GA found
fval = inf;
for i = 1:length(ga_state)
    % go through all runs and find global maxima
    [cur_opt, cur_ind] = min(ga_state{i}.score);
    if cur_opt < fval
        fval = cur_opt;
        histogram = ga_state{i}.population(cur_ind, :).';
    end
end

% get histogram which is really used by my Histogram function
[~, used_hist] = Histogram_Importance_Sampling(1:10, arch.us.num_samples, arch.us.hist.boundaries, 'given', histogram);
histogram = used_hist;
optim_nmse = fval;

fprintf('Best histogram [');
fprintf('%d, ', histogram);
fprintf('] with NMSE = %f\n', fval);
