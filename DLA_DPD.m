function [tx_signal_pa, coefs, res] = DLA_DPD(tx_signal, coefs, pars)

% calculate DPD output
tx_signal_norm = tx_signal / pars.pa_model.in_scaler;
U = pars.dpd.Matrix_Func(tx_signal_norm, pars.dpd.P, pars.dpd.M);
dpd_out_norm = U * coefs;
dpd_out = dpd_out_norm * pars.pa_model.in_scaler;

% calculate PA model output
%tx_signal_pa = pars.PA_Model_func(dpd_out, pars.pa.coefs, pars.pa.P, pars.pa.M);
tx_signal_pa = pars.pa_model.GetOutput(dpd_out);

if isempty(tx_signal_pa)
    res.nmse = 10;
    res.cond = NaN;
    return;
end

tx_signal_ph = exp(-1j*(pars.pa.ph)) * tx_signal;   % phase needs to be approx. synchronised
tx_signal_pan = awgn(tx_signal_pa, pars.fb.SNR, db(0.5,'power')); % add noise to feedback signal

% scale signals for DPD calculation
tx_signal_pan_norm = tx_signal_pan / pars.pa_model.out_scaler / pars.pa.add_gain;
tx_signal_ph = tx_signal_ph / pars.pa_model.in_scaler;

% gain = lscov(tx_signal_ph, tx_signal_pa);
% sprintf('abs(gain) = %.8f\nangle(gain) = %.8f', abs(gain), angle(gain))
% % in Prepare_PA_Model
% error('konec');

% undersampling if specified
if isfield(pars.dpd, 'us') && pars.dpd.us.enable
    if strcmp(pars.dpd.us.method, 'histogram') || strcmp(pars.dpd.us.method, 'hist')
        us = Histogram_Importance_Sampling(abs(dpd_out_norm), pars.dpd.us.num_samples, ...
            pars.dpd.us.hist.boundaries, pars.dpd.us.hist.method, pars.dpd.us.hist.reference);

        % check if histogram method provided non-zero length of the signal
        % if so, do not update coefficients and return high NMSE
        if length(dpd_out(us)) < pars.dpd.us.num_samples - 5
            res.nmse = 10;
            res.cond = NaN;
            return;
        end
    elseif strcmp(pars.dpd.us.method, '1N')
        us_fact = ceil(length(dpd_out)/pars.dpd.us.num_samples);
        us = 1:us_fact:length(dpd_out);
    elseif strcmp(pars.dpd.us.method, 'block')
        start = round((length(dpd_out)-pars.dpd.us.num_samples)/2);
        us = start:start+pars.dpd.us.num_samples-1;
    elseif strcmp(pars.dpd.us.method, 'qr')
        [~, us] = QR_Reduce_Matrix([pars.dpd.Matrix_Func(tx_signal_ph, pars.dpd.P, pars.dpd.M), tx_signal_pan_norm - tx_signal_ph].', ...
            pars.dpd.us.num_samples);
    elseif strcmp(pars.dpd.us.method, 'gs')
        us = Gradient_Sampling(pars.dpd.Matrix_Func(tx_signal_ph, pars.dpd.P, pars.dpd.M), ...
            tx_signal_pan_norm - tx_signal_ph, pars.dpd.us.num_samples);
    elseif strcmp(pars.dpd.us.method, 'gs_qr')
        us = Gradient_Sampling(pars.dpd.Matrix_Func(tx_signal_ph, pars.dpd.P, pars.dpd.M), ...
            tx_signal_pan_norm - tx_signal_ph, pars.dpd.us.num_samples+100);
        [~, us2] = QR_Reduce_Matrix(pars.dpd.US_Matrix_Func(tx_signal_ph, pars.dpd.P, pars.dpd.M, us).', ...
            pars.dpd.us.num_samples);
        us3 = zeros(sum(us),1);
        us3(us2) = 1;
        us3 = us3 > 0;
        idx = find(us > 0);
        us(idx(~us3)) = 0;
%         fprintf('Number of GS QR samples = %d\n', sum(us));

        if isfield(pars.dpd.us, 'hist_training') && pars.dpd.us.hist_training
            % save samples into results and plot histogram here
            res.sel_samples = dpd_out_norm(us);

            if isfield(pars.dpd.us, 'is_hist_plot') && pars.dpd.us.is_hist_plot
                % plot histograms
                figure(21);
                histogram(abs(dpd_out_norm), pars.dpd.us.hist.boundaries);
                xlim([0 1]);

                figure(22);
                histogram(abs(dpd_out_norm(us)), pars.dpd.us.hist.boundaries);
                xlim([0 1]);
                histcounts(abs(dpd_out_norm(us)), pars.dpd.us.hist.boundaries)
            end
        end
    else
        error('Undersampling method %s is not supported', pars.dpd.us.method);
    end
%     dla_samples = length(dpd_out(us))
else
    us = 1:length(dpd_out);
end

% use Gauss-Newton method to update coefficients
U = pars.dpd.US_Matrix_Func(tx_signal_ph, pars.dpd.P, pars.dpd.M, us);
coefs = coefs - pars.dpd.dla.mu*((U'*U)\(U'*(tx_signal_pan_norm(us)-tx_signal_ph(us))));
% coefs = coefs - pars.dpd.dla.mu*lscov(U,tx_signal_pan_norm(us)-tx_signal_ph(us));

% if strcmp(pars.dpd.us.method, 'histogram') || strcmp(pars.dpd.us.method, 'hist')
    % plot AM/AM characteristisc
    figure(31);
    subplot(2,1,1);
    plot(abs(tx_signal)/pars.pa_model.in_scaler, abs(tx_signal_pa), '.',  'DisplayName','Linearised');
    hold on;
    plot(abs(tx_signal)/pars.pa_model.in_scaler, abs(dpd_out)/pars.pa_model.in_scaler*pars.pa_model.out_scaler, '.',  'DisplayName','DPD out');
    plot(abs(dpd_out)/pars.pa_model.in_scaler, abs(tx_signal_pa), '.', 'DisplayName','PA');
    plot(abs(dpd_out(us))/pars.pa_model.in_scaler, abs(tx_signal_pa(us)), '+', 'DisplayName', 'Selected');
    hold off;
    title('AM/AM characteristics');
    legend('show', 'Location','southeast');
    subplot(2,1,2);
    ph_diff = angle(tx_signal_pa) - angle(tx_signal); ph_diff(ph_diff > pi) = ph_diff(ph_diff > pi) - 2*pi; ph_diff(ph_diff < -pi) = ph_diff(ph_diff < -pi) + 2*pi;
    plot(abs(tx_signal), ph_diff, '.', 'DisplayName','Linearised');
    hold on;
    ph_diff = angle(dpd_out) - angle(tx_signal); ph_diff(ph_diff > pi) = ph_diff(ph_diff > pi) - 2*pi; ph_diff(ph_diff < -pi) = ph_diff(ph_diff < -pi) + 2*pi;
    plot(abs(tx_signal), ph_diff, '.', 'DisplayName','DPD out');
    ph_diff = angle(tx_signal_pa) - angle(dpd_out); ph_diff(ph_diff > pi) = ph_diff(ph_diff > pi) - 2*pi; ph_diff(ph_diff < -pi) = ph_diff(ph_diff < -pi) + 2*pi;
    plot(abs(dpd_out), ph_diff, '.', 'DisplayName','PA');
    title(sprintf('%s, num samples = %d', pars.dpd.us.method, length(tx_signal_ph(us))));
    hold off;
% end

% calculate condition number of the inversed matrix
res.cond = cond(U'*U);

% calculate NMSE
tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
res.nmse = nmse(tx_signal, tx_signal_pan);

end
