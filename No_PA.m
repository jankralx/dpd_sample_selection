function [tx_signal_pa, coefs, res] = No_PA(tx_signal, coefs, pars)

% no PA, input goes directly to the output
tx_signal_pa = tx_signal;

% calculate NMSE
res.nmse = inf;
res.cond = 0;

end
