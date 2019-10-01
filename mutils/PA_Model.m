classdef PA_Model
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pa_descript_fname = 'pa_models/pa_descript.mat';
        pa_name;
        fc = 1e9;
        fs = 10e6;
        modulation = 'QAM16';
        coefs;
        cur_coefs;
        
        model_name;
        pa_coefs;
        in_scaler;
        out_scaler;
        mpars;
        pa_meas;
        max_pa_peak = -30;          % in dBm
    end
    
    methods
        function obj = PA_Model(opt)

%             if nargin >= 1 && isfield(opt,'pa_name')
%                 obj.pa_name = opt.pa_name;
%             end
%             if nargin >= 1 && isfield(opt,'fc')
%                 obj.fc = opt.fc;
%             end
%             if nargin >= 1 && isfield(opt,'fs')
%                 obj.fs = opt.fs;
%             end
            
            % temporary solution
            if nargin >= 1 && isfield(opt,'coef_fname')
                load(opt.coef_fname);
                obj.model_name = model_name;
                obj.pa_coefs = pa_coefs;
                obj.in_scaler = in_scaler;
                obj.out_scaler = out_scaler;
                obj.mpars = mpars;
            end
            
            if nargin >= 1 && isfield(opt,'pa_meas')
                obj.pa_meas = opt.pa_meas;
            end
            
            if exist('pa_exp_region', 'var')
                obj.pa_meas.pa_exp_region = pa_exp_region;
            end
            
            if nargin >= 1 && isfield(opt,'max_pa_peak')
                obj.max_pa_peak = opt.max_pa_peak;
            end
                        
            % based on the given parameters find the PA measurement in the pa_description file
            %load(obj.pa_descript_fname);
            
%             for i = 1:length(pa_model_params)
%             end        
        end
        
        
        function output = GetOutput(obj, input)
            
            % calculate model output based on used model
            switch obj.model_name
                case 'GMP'
                    % scale the input
                    input = input / obj.in_scaler;
                    % calculate model output
                    output = GMP_Output(input, obj.pa_coefs, obj.mpars.Ka, ...
                        obj.mpars.La, obj.mpars.Kb, obj.mpars.Lb, obj.mpars.Mb, ...
                        obj.mpars.Kc, obj.mpars.Lc, obj.mpars.Mc);
                    % scale output
                    output = output * obj.out_scaler;
                case 'DDR2'
                    % scale the input
                    input = input / obj.in_scaler;
                    % calculate model output
                    output = DDR2_Output(input, obj.mpars.P, obj.mpars.M, obj.pa_coefs);
                    % scale output
                    output = output * obj.out_scaler;
                case 'MP'
                    % scale the input
                    input = input / obj.in_scaler;
                    % calculate model output
                    output = MP_Output(input, obj.mpars.P, obj.mpars.M, obj.pa_coefs);
                    % scale output
                    output = output * obj.out_scaler;
                case 'RF_WebLab'
                    if ~isvector(input) || nnz(isnan(input)) || nnz(abs(input) == Inf)
                        output = [];
                        return;
                    end
                    
                    % limit input to -11 dBm
                    R = 50;
                    dblim = -11;
                    vlim = sqrt(10^(dblim/10)/1000*R);
                    over = abs(input) > vlim;
                    if nnz(over)
                        % clip values greater than -11dBm, preserve phase
                        input(over) = vlim * exp(1j*angle(input(over)));
                        warning('Clipping values greater than %.1f dBm', dblim);
                    end
                    
                    try
                        output = RF_WebLab_GetOutput(input);
                    catch
                        output = [];
                    end
                case 'SMU200A_FSVR'
                    if ~isvector(input) || nnz(isnan(input)) || nnz(abs(input) == Inf)
                        output = [];
                        return;
                    end

                    % limit input to specified dBm value
                    R = 50;
                    vlim = sqrt(10^(obj.max_pa_peak/10)/1000*R);
                    over = abs(input) > vlim;
                    if nnz(over)
                        % clip values greater than defined limit, preserve phase
                        input(over) = vlim * exp(1j*angle(input(over)));
                        warning('Clipping values greater than %.1f dBm', obj.max_pa_peak);
                    end
                    
                    output = obj.pa_meas.GetOutput(input);
                        
                otherwise
                    error('Model name %s not supported', obj.model_name);
            end
            
        end
    end
end

