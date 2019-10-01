function [tx_signal_pa, coefs, res] = No_DPD(tx_signal, coefs, pars)

% add PA backoff
tx_signal = tx_signal / 10^(pars.pa.backoff/20);

% transmit the signal through the PA model
tx_signal_pa = pars.pa_model.GetOutput(tx_signal);

% calculate NMSE
tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
res.nmse = nmse(tx_signal, tx_signal_pan);
res.cond = 0;

end
