function criteria=HistogramFitness(hist, pars, arch)

% if architecture has us (undersampling) field, then copy it to
% local parameters
if isfield(arch,'us')
    pars.dpd.us = arch.us;
else
    error('Cannot optimise histogram when undersampling is disabled');
end

pars.dpd.us.hist.reference = hist;

if pars.histopt.num_iter < 10
    error('number of itterationts too low');
end

criteria = [];
% make multiple runs to include probability of convergence
for run = 1:pars.histopt.num_runs
    % make few iterataions to evaluate criterion function
    for iter = 1:pars.histopt.num_iter
        % update mu
        pars.dpd.dla.mu = pars.dpd.dla.mu0 / pars.dpd.dla.gamma.^(...
        (iter - pars.dpd.dla.change_mu_after)*...
        (iter > pars.dpd.dla.change_mu_after));

        % create tx signal for transmission
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
        % adjust transmitt power
        pep_out = pars.pa.pep;          % output peak power in dBm
        tx_sig = tx_signal / max(abs(tx_signal)) * (10^((pep_out-30)/20)*sqrt(50));

        % calculate iteration of DPD architecture
        coefs = arch.coefs(iter,:).';
        [sig_pa, coefs, res] = arch.func(tx_sig, coefs, pars);

        % calculate ACPRs and power in the main channel
        [acpr_val, mchan_pwr] = ACPR(sig_pa, ...
            pars.mod.oversample_fact, pars.res.acrp.channels, false);

        % save results locally
        arch.acpr(iter,:,:) = acpr_val;
        arch.mchan_pwr(iter) = mchan_pwr;

        arch.res(iter) = res;
        arch.coefs(iter+1,:) = coefs.';
    end

    % calculate output criterion function as average of a certain number of
    % last iterations
    nmserr = [arch.res(:).nmse];

    criteria(run) = Avg_dB(nmserr(pars.histopt.avg_after:pars.histopt.num_iter).', 1, 10);
    if isnan(criteria)
        criteria(run) = 10;
    end

    % figure(1);
    % plot(nmserr);

end
criteria = Avg_dB(criteria.', 1, 10);





