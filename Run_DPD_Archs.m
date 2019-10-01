function arch = Run_DPD_Archs(pars, arch)

freq_axis = linspace(-0.5, 0.5, pars.burst_len + 1);
freq_axis = freq_axis(2:end);
arch_num = size(arch,1);

% initialize architecture fields
if pars.is_fft
    for i = 1:arch_num
        arch(i).fftPA_sig = zeros(pars.burst_len,1);
        arch(i).fftPA_sig_avg = zeros(pars.burst_len,1);
        arch(i).fftPA_sig_max = zeros(pars.burst_len,1);
    end
end

for iter = 1:pars.num_iter
%     fprintf('Iter = %i\n',iter);
    % set mu based on the iteration
    pars.dpd.dla.mu = pars.dpd.dla.mu0 / pars.dpd.dla.gamma.^(...
        (iter - pars.dpd.dla.change_mu_after)*...
        (iter > pars.dpd.dla.change_mu_after));
    if pars.dpd.dla.mu < pars.dpd.dla.mu_min
        pars.dpd.dla.mu = pars.dpd.dla.mu_min;
    end
    % create tx signal for transmission
    if strcmp(pars.mod.type, 'QAM')
        binin = rand(pars.mod.num_bits, 1) > 0.5;
        tx_signal = QAM_Modulator(binin, pars.mod);
    elseif strcmp(pars.mod.type, 'OFDM')
        bits_per_frame = log2(pars.mod.imod.M)*pars.mod.rb_size*pars.mod.num_rb;
        tx_bin = randi([0, 1], pars.mod.num_frames*bits_per_frame, 1);           % generate random binary data
        [tx_signal, tx_symbs] = OFDM_Modulator(tx_bin, pars.mod);
    else
        tx_signal = FBMC_modulator(pars.mod.num_subchannels, ...
            pars.mod.oversample_fact, pars.mod.num_frames);
    end

    % adjust transmitt power
    pep_out = pars.pa.pep;          % output peak power in dBm
    tx_signal = tx_signal / max(abs(tx_signal)) * (10^((pep_out-30)/20)*sqrt(50));
    
    for i = 1:arch_num
        tx_sig = tx_signal;
        locpars = pars;
        
        % if architecture has us (undersampling) field, then copy it to
        % current parameters
        if isfield(arch(i),'us')
            locpars.dpd.us = arch(i).us;
        end
        
        coefs = arch(i).coefs(iter,:).';
        [sig_pa, coefs, res] = arch(i).func(tx_sig, coefs, locpars);
        
        try
            % calculate ACPRs and power in the main channel
            [acpr_val, mchan_pwr] = ACPR(sig_pa, ...
                pars.mod.oversample_fact, pars.res.acrp.channels, false);
        catch
            acpr_val = NaN * ones(2);
            mchan_pwr = NaN;
        end

%         if i == length(arch_names)   % backoff only for PA without DPD
%             backoff = 10^((mchan_pwr-arch(1).mchan_pwr(iter))/20);
%             [sig_pa, err, coefs] = arch(i).func(tx_signal/backoff, coefs, pars);
%         
%             % calculate ACPRs and power in the main channel
%             [acpr_val, mchan_pwr] = ACPR(sig_pa, ...
%                 pars.mod.oversample_fact, pars.res.acrp.channels, false);
%         end

        % calculate and save EVM
        if strcmp(pars.mod.type, 'OFDM')
            try
                [~, rx_symbs] = OFDM_Demodulator(sig_pa, pars.mod, tx_symbs);   % demodulation
                evm = EVM(rx_symbs, tx_symbs);                                  % calculate EVM
            catch
                evm = inf;
            end
        else
            evm = 0;
        end

        % save parameters
        arch(i).evm(iter) = evm;
        arch(i).acpr(iter,:,:) = acpr_val;
        arch(i).mchan_pwr(iter) = mchan_pwr;

        arch(i).res(iter) = res;
        arch(i).coefs(iter+1,:) = coefs.';
        
        if pars.is_fft
            if pars.is_plot || iter > pars.res.avg_fft_after_iter
                try
                    arch(i).fftPA_sig = fft(sig_pa);       % calculate spectrum if it is needed
                catch
                    arch(i).fftPA_sig = ones(size(tx_signal));
                end
            end

            if iter > pars.res.avg_fft_after_iter
                arch(i).fftPA_sig_avg = arch(i).fftPA_sig_avg + abs(arch(i).fftPA_sig);
                arch(i).fftPA_sig_max = max([arch(i).fftPA_sig_max, abs(arch(i).fftPA_sig)],[],2);
            end
        end
    end
    
    % all architectures have been evaluated, plot the results now
    if pars.is_plot
        % plot amplitude frequency spectra
        figure(1);
        for i = 1:size(arch,1)
            if i == 2; hold on; end
            try
                plot(freq_axis, db(fftshift(arch(i).fftPA_sig)), 'DisplayName',arch(i).name);
            catch
            end
        end
        hold off;
        legend('show');
        title('Magnitude spectra');
        xlabel('Frequency');
        ylabel('Magnitude (dB)');
        axis([-0.5 0.5, -100 60]);

        % NMSE evolution through iterations
        figure(2);
        ylimits = [inf -inf];
        for i = 1:size(arch,1)
            if i == 2; hold on; end
            plot([arch(i).res(1:iter).nmse], '-x', 'DisplayName',arch(i).name);
            
            miny = min([arch(i).res(1:iter).nmse]);
            maxy = max([arch(i).res(1:iter).nmse]);
            if ylimits(1) > miny
                ylimits(1) = round(miny-10,-1);
            end
            if ylimits(2) < maxy
                ylimits(2) = round(maxy+10,-1);
            end
        end
        hold off;
        legend('show');
        title('Normalised mean square error');
        xlabel('Iteration cycle');
        ylabel('NMSE (dB)');
        axis([0 pars.num_iter, ylimits(1), ylimits(2)]);
        
        % ACPR evolution through iterations
        figure(3);
        for i = 1:size(arch,1)
            if i == 2
                subplot(2,1,1);
                hold on;
                subplot(2,1,2);
                hold on;
            end
            subplot(2,1,1);
            plot(Avg_dB(arch(i).acpr(1:iter,1,:),3,10), '-x', 'DisplayName',arch(i).name);

            subplot(2,1,2);
            plot(Avg_dB(arch(i).acpr(1:iter,2,:),3,10), '-x', 'DisplayName',arch(i).name);
        end
        subplot(2,1,1);
        hold off;
        legend('show');
        title('ACPR in 1-st adjacent channel');
        xlabel('Iteration cycle');
        ylabel('ACPR (dB)');
        axis([0 pars.num_iter, -70 0]);
        
        subplot(2,1,2);
        hold off;
        legend('show');
        title('ACPR in 2-nd adjacent channel');
        xlabel('Iteration cycle');
        ylabel('ACPR (dB)');
        axis([0 pars.num_iter, -70 0]);
        
        % Main channel power
        figure(4);
        ylimits = [inf -inf];
        for i = 1:size(arch,1)
            if i == 2; hold on; end
            plot(arch(i).mchan_pwr(1:iter), '-x', 'DisplayName',arch(i).name);
            
            miny = min(arch(i).mchan_pwr(1:iter));
            maxy = max(arch(i).mchan_pwr(1:iter));
            if ylimits(1) > miny
                ylimits(1) = round(miny-10,-1);
            end
            if ylimits(2) < maxy
                ylimits(2) = round(maxy+10,-1);
            end
            
            fprintf('Mean output power in main channel = %.3f dBW\n', mean(arch(i).mchan_pwr(1:iter)));
        end
        hold off;
        legend('show');
        title('Main channel power');
        xlabel('Iteration cycle');
        ylabel('P_{mch} (dBW)');
        axis([0 pars.num_iter, ylimits(1), ylimits(2)]);

        % EVM
        figure(5);
        ylimits = [inf -inf];
        for i = 1:size(arch,1)
            if i == 2; hold on; end
            plot(arch(i).evm(1:iter)*100, '-x', 'DisplayName',arch(i).name);
            
            miny = min(arch(i).evm(1:iter)*100);
            maxy = max(arch(i).evm(1:iter)*100);
            if ylimits(1) > miny
                ylimits(1) = round(miny-10,-1);
            end
            if ylimits(2) < maxy
                ylimits(2) = round(maxy+10,-1);
            end
        end
        hold off;
        legend('show');
        title('EVM');
        xlabel('Iteration cycle');
        ylabel('EVM (%)');
        axis([0 pars.num_iter, ylimits(1), ylimits(2)]);
        
        % Condition Number of Matrix for Inversion
        figure(51);
        for i = 1:size(arch,1)
            if i == 2; hold on; end
            semilogy([arch(i).res(1:iter).cond], '-x', 'DisplayName',arch(i).name);
        end
        hold off;
        legend('show');
        title('Condition Number of Inversed Matrix');
        xlabel('Iteration cycle');
        ylabel('Condition number (-)');
        axis([0 pars.num_iter, -inf inf]);
    end
    drawnow;
end
%%
if pars.num_iter > pars.res.avg_fft_after_iter
    % average saved spectra
    for i = 1:size(arch,1)
        arch(i).fftPA_sig_avg = arch(i).fftPA_sig_avg / (pars.num_iter - pars.res.avg_fft_after_iter);
    end

    if pars.is_plot
        % plot result spectra
        figure(1);
        for i = 1:size(arch,1)
            if i == 2; hold on; end
            plot(freq_axis, db(fftshift(arch(i).fftPA_sig_avg)), 'DisplayName',arch(i).name);
        end
        hold off;
        legend('show');
        title('Magnitude spectra');
        xlabel('Frequency');
        ylabel('Magnitude (dB)');
        axis([-inf inf, -30 60]);

    end
end
