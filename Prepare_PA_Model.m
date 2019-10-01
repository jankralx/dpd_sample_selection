% PA coefficients selection
switch pars.pa.sel

    case 18
        % Hon Driver
        pa_model_opt.coef_fname = 'amplifiers/HonDri_OFDM_GMP_coefs.mat';
        pars.pa_model = PA_Model(pa_model_opt);
        pars.pa.pep = 4;
        pars.pa.ph = 0;
        pars.pa.add_gain = 1.3;
        
end
