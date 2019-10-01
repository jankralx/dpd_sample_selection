addpath('mutils');
addpath('histogram_optim');

mkdir('hist_results');

% clear all;

global Fs;
Fs = 30e6;

%--------------------------------------------------------------------------
% Simulation parameters
%--------------------------------------------------------------------------
pars.is_test = 0;
pars.is_hist_training = 0;
pars.is_no_dpd_calc = 0;

% PA model selection
pars.pa.sel = 18;

% set to true if plots are required
pars.is_plot = 0;

% set to true if FFTs are to be saved in results
pars.is_fft = pars.is_plot;

% FBMC Modulator
pars.fbmc.type = 'FBMC';
pars.fbmc.num_subchannels = 1024; % num of SUBCHANNELS
pars.fbmc.oversample_fact = 12; % 12 % oversampling
pars.fbmc.num_frames = 26; % 18

% QAM Modulator
pars.qam.type = 'QAM';
pars.qam.num_bits = 1024*8;
pars.qam.M = 16;
pars.qam.oversample_fact = 5; % 12 % oversampling
pars.qam.filt = 'sqrt';

% f-OFDM modulator
pars.ofdm.type = 'OFDM';
pars.ofdm.num_frames = 12;
pars.ofdm.imod.M = 64;
pars.ofdm.num_rb = 68;
pars.ofdm.rb_size = 12;
pars.ofdm.fft_len = 4096;
pars.ofdm.filt.L = pars.ofdm.fft_len/2 + 1;
pars.ofdm.filt.ena = true;      % enable filter for f-OFDM
pars.ofdm.oversample_fact = 5;   % because of ACPR calculation



% Modulator selection
pars.mod = pars.ofdm;


% create tx signal -  we need dimensions of frequency axis
if strcmp(pars.mod.type, 'QAM')
    binin = rand(pars.mod.num_bits, 1) > 0.5;
    tx_signal = QAM_Modulator(binin, pars.mod);
elseif strcmp(pars.mod.type, 'OFDM')
    bits_per_frame = log2(pars.mod.imod.M)*pars.mod.rb_size*pars.mod.num_rb;
    tx_bin = randi([0, 1], pars.mod.num_frames*bits_per_frame, 1);           % generate random binary data
    [tx_signal, ~] = OFDM_Modulator(tx_bin, pars.mod);
else
    tx_signal = FBMC_modulator(pars.mod.num_subchannels, ...
        pars.mod.oversample_fact, pars.mod.num_frames);
end
pars.burst_len = length(tx_signal);


% noise in feedback
pars.fb.SNR = inf;

% number of repetitions for convergence probability calculation
pars.num_prob_rep = 1000;
% number of repetitions in Histogram Fitness function
pars.histopt.num_runs = 10;
% number of generations used by GA
pars.histopt.generations = 20;
% population size used by GA
pars.histopt.popsize = 100;

% histogram optimization parameters
pars.histopt.num_iter = 80;
pars.histopt.avg_after = 60;

% Number of iterations
pars.num_iter = 80;
pars.res.avg_fft_after_iter = 60;        % average FFT results after specified number of iterations

pars.res.acrp.channels = [1/2+0.1, 3/2+0.1,...
                          3/2+0.2, 5/2+0.2];

%
% DPD parameters
%

pars.dpd.Matrix_Func = @(tx_sig, P, M)(DDR2_Matrix(tx_sig, P, M));
pars.dpd.US_Matrix_Func = @(tx_sig, P, M, us)(DDR2_UndersampledMatrix(tx_sig, P, M, us));
pars.dpd.P = 7;
pars.dpd.M = 1;

% pars.dpd.Matrix_Func = @(tx_sig, P, M)(DVR_Matrix(tx_sig, P, M));
% pars.dpd.US_Matrix_Func = @(tx_sig, P, M, us)(DVR_UndersampledMatrix(tx_sig, P, M, 3, us));
%
% % pars.dpd.P = [0 0.01 0.065 0.17 0.3 0.6 1];
% % pars.dpd.M = 5;
% pars.dpd.P = [0 0.065 0.3 1];
% pars.dpd.M = 2;


% DLA parameters
pars.dpd.dla.mu0 = 0.05;
pars.dpd.dla.change_mu_after = 15;
pars.dpd.dla.mu_min = 0.01;

% gamma is calculated such it makes mu to be mu_min after pars.res.avg_fft_after_iter
% iterations
pars.dpd.dla.gamma = (pars.dpd.dla.mu0/pars.dpd.dla.mu_min).^...
    (1/(pars.res.avg_fft_after_iter-pars.dpd.dla.change_mu_after));

if pars.is_test
    % number of repetitions for convergence probability calculation
    pars.num_prob_rep = 1;
    % number of repetitions in Histogram Fitness function
    pars.histopt.num_runs = 1;
    % number of generations used by GA
    pars.histopt.generations = 5;
    % population size used by GA
    pars.histopt.popsize = 4;
end
%------------------------------------------------------------------------------
% Define Architectures
%------------------------------------------------------------------------------

pars.num_arch = 7;

% ___default architecture parameters___

% initialize coefficients such first coefficient is 1 and others are 0
dpd_num_coefs = size(pars.dpd.Matrix_Func((1:10).', pars.dpd.P, pars.dpd.M), 2);
def_coefs = zeros(dpd_num_coefs,1);
def_coefs(1) = 0.5;

def.name = '';
def.func = @(tx_signal, coefs, dpd_par)(No_DPD(tx_signal, coefs, dpd_par));
def.coefs = zeros(pars.num_iter+1,length(def_coefs));
def.coefs(1,:) = def_coefs;
def.fftPA_sig_max = [];
def.fftPA_sig_avg = 0;
def.evm = zeros(pars.num_iter,1);
def.mchan_pwr = zeros(pars.num_iter,1);
def.acpr = zeros(pars.num_iter,2,2);

def.us.enable = true;
def.us.method = 'block';                 % {'block', '1N', 'histogram'}
def.us.num_samples = 0;                  % target number of points after undersampling

% create architectures without DPD and without PA to evaluate enhancements
% caused by the DPD
arch_no_pa = def;
arch_no_pa.func = @(tx_signal, coefs, dpd_par)(No_PA(tx_signal, coefs, dpd_par));
arch_no_dpd = def;

% create parameter array from the default values
arch = repmat(def,pars.num_arch,1);
clear def;   % clear default values as they will not be needed anymore

% ___set individual architectures___

i = 1;

% % classical ILA with block undersampling (N continous samples taken)
% arch(i).name = 'ILA';
% arch(i).func = @(tx_signal, coefs, pars)(ILA_DPD(tx_signal, coefs, pars));
% i = i + 1;

% DLA with block undersampling (N continous samples taken)
arch(i).name = 'Block';
arch(i).func = @(tx_signal, coefs, pars)(DLA_DPD(tx_signal, coefs, pars));
i = i + 1;

% DLA with 1/N undersampling
arch(i).name = '1/N';
arch(i).func = @(tx_signal, coefs, pars)(DLA_DPD(tx_signal, coefs, pars));
arch(i).us.method = '1N';
i = i + 1;


% DLA with even histogram
arch(i).name = 'Hist-Even';
arch(i).func = @(tx_signal, coefs, pars)(DLA_DPD(tx_signal, coefs, pars));
arch(i).us.method = 'hist';
arch(i).us.hist.method = 'even';
arch(i).us.hist.boundaries = [0:0.1:0.9, 2].';
arch(i).us.hist.reference = ones(length(arch(i).us.hist.boundaries)-1,1);
i = i + 1;

% DLA with given histogram - optimised by GA
arch(i).name = 'Hist-GA';
arch(i).func = @(tx_signal, coefs, pars)(DLA_DPD(tx_signal, coefs, pars));
arch(i).us.method = 'hist';
arch(i).us.hist.method = 'given';
arch(i).us.hist.boundaries = [0:0.1:0.9, 2].';
arch(i).us.hist.reference = ones(length(arch(i).us.hist.boundaries)-1,1);
arch(i).us.hist.optim_nmse = Inf;
i = i + 1;

% DLA with Tom Gotthans QR decomposition
arch(i).name = 'QR';
arch(i).func = @(tx_signal, coefs, pars)(DLA_DPD(tx_signal, coefs, pars));
arch(i).us.method = 'qr';
i = i + 1;

% DLA with Gradient Sampling
arch(i).name = 'GS';
arch(i).func = @(tx_signal, coefs, pars)(DLA_DPD(tx_signal, coefs, pars));
arch(i).us.method = 'gs';
i = i + 1;

% DLA with Gradient Sampling + QR
arch(i).name = 'GS-QR';
arch(i).func = @(tx_signal, coefs, pars)(DLA_DPD(tx_signal, coefs, pars));
arch(i).us.method = 'gs_qr';
i = i + 1;



% % Without DPD (supressed for final runs - not needed as there is no
% % undersampling involved)
% arch(i).name = 'No DPD';
% arch(i).func = @(tx_signal, coefs, pars)(No_DPD(tx_signal, coefs, pars));
% i = i + 1;

%--------------------------------------------------------------------------
% Run architectures
%--------------------------------------------------------------------------
if ~pars.is_test
    p = gcp;
    warning('off', 'all');
    pctRunOnAll warning('off', 'all');
end

Prepare_PA_Model;       % selects PA model

if pars.is_test
    num_samples_simrange = 1000;
else
    num_samples_simrange = [17 18 19 22 24 28 32 50 100 200 400 700 1000 2000 20000 50000];
end

% prepare cell array for results
arch_res = cell(length(num_samples_simrange), pars.num_prob_rep);
arch_hist = cell(length(num_samples_simrange),1);

if pars.is_hist_training
    for num_sampl_ind = 1:length(num_samples_simrange)
        fprintf('Runing histogram optimisation %d/%d for %i samples\n', num_sampl_ind, length(num_samples_simrange), num_samples_simrange(num_sampl_ind));
        tic
        % optimize histogram for all architectures using histogram
        for i = 1:size(arch,1)
            arch(i).us.num_samples = num_samples_simrange(num_sampl_ind);
            if (strcmp(arch(i).us.method, 'histogram') || strcmp(arch(i).us.method, 'hist')) && ...
                    strcmp(arch(i).us.hist.method, 'given')
                % optimize histogram for this architecture
                [arch(i).us.hist.reference, arch(i).us.hist.optim_nmse] = Optimize_Histogram(pars, arch(i));
            end
        end
        toc

        arch_hist{num_sampl_ind} = arch;
        save(sprintf('hist_results/results_01_hist_%d.mat',num_sampl_ind), 'arch_hist', 'arch', 'pars', 'num_samples_simrange', '-v7.3');
    end

    save('results_01_hist.mat', 'arch_hist', 'arch', 'pars', 'num_samples_simrange', '-v7.3');
else
    % load architectures with already optimised histograms
    load('results_01_hist.mat', 'arch_hist');
end

if ~pars.is_no_dpd_calc
    for num_sampl_ind = 1:length(num_samples_simrange)
        for i = 1:size(arch,1)
            arch(i).us.num_samples = num_samples_simrange(num_sampl_ind);
        end
        % find architecture histogram which should be used for evaluation of
        % num_samples_simrange(num_sampl_ind)
        opt_found = 0;
        for i = 1:length(arch_hist)
            if arch_hist{i}(1).us.num_samples == num_samples_simrange(num_sampl_ind) || (num_samples_simrange(num_sampl_ind) >= 2000 && arch_hist{i}(1).us.num_samples == 2000)
                arch(4).us.hist.reference = arch_hist{i}(4).us.hist.reference;
                opt_found = 1;
            end
        end
        if ~opt_found
            error('Could not find optimized histogram for %i samples', ...
                num_samples_simrange(num_sampl_ind));
        end

        ares = cell(1, pars.num_prob_rep);
        % run multiple times to get some statistics
        if pars.is_test
            fprintf(' Runing evaluation %d/%d for %i samples\n', ...
                num_sampl_ind, length(num_samples_simrange), num_samples_simrange(num_sampl_ind));
            for rep = 1:pars.num_prob_rep
                ares{rep} = Run_DPD_Archs(pars, arch);
            end
        else
            fprintf(' Runing evaluation %d/%d for %i samples in parallel\n', ...
                num_sampl_ind, length(num_samples_simrange), num_samples_simrange(num_sampl_ind));
            parfor rep = 1:pars.num_prob_rep
                ares{rep} = Run_DPD_Archs(pars, arch);
            end
        end

        % save results from parallel runs into output results
        for rep = 1:pars.num_prob_rep
            arch_res{num_sampl_ind, rep} = ares{rep};
        end
    end

    % save results
    save('results_01.mat', 'arch_res', 'arch_hist', 'arch', 'pars', 'num_samples_simrange', '-v7.3');
else
    % calculate results for selected PA model without DPD
    pars_one_iter = pars;
    pars_one_iter.pa.backoff = 0.2;                  % in dB
    ares_no_dpd = Run_DPD_Archs(pars_one_iter, arch_no_dpd);
    ares_no_pa = Run_DPD_Archs(pars_one_iter, arch_no_pa);

    % save results
    save('results_01_nodpd.mat', 'ares_no_dpd', 'ares_no_pa', 'pars', '-v7.3');
end
