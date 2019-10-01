function [ params ] = RegisterParam( params, par_name, def_val )
%REGISTERPARAM Register parameter of a function.
%   The function check if given name is a field of the given structure and
%   if it is not, it sets it to the given default value.
%
%   If no def_val is given then empty struct is created for given parameter
%   name.

if ~isfield(params, par_name)
    if nargin >= 3
        params.(par_name) = def_val;
    else
        params.(par_name) = struct();
    end
end

end

