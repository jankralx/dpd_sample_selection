# Feedback Sample Selection Methods Allowing Lightweight Digital
# Predistorter Adaptation

The Matlab source codes are provided to supplement our research paper
"Feedback Sample Selection Methods Allowing Lightweight Digital
Predistorter Adaptation". They allows to reproduce simulation results
in our paper.

## Requirements

We have run simulations on Ubuntu OS, Matlab 2018a, but the
simulations should be OS independent and all Matlab versions > 2018 should
be compatible. Please note that for higher simulation executation
Distributed Computation Toolbox is required, however, all simulations
can be done without parfor loops and hence should require no
toolboxes.

## Simulation Execution

The source codes are split into two parts: 1) calculation of
results, 2) generating plots.

1) run *RUN_ANALYSIS_04.m* to generate simulation results
2) run *PLOT_RESULTS_04.m* to plot the simulation results

Please note that simulations are quite time demanding and can execute
more than 24 hours (depending on the performance of HW).

## Extended Simulations

By default, Matlab script *RUN_ANALYSIS_04.m* loads already optimised
histograms from *results_01_hist.mat*. You can turn on histogram
optimisation by changing line 15 of *RUN_ANALYSIS_04.m* to
"pars.is_hist_training = 1;".


## Please Cite Our Paper
