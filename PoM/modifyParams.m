function [param_vals_scaled, modParam_baseline, modParam_vals]  = modifyParams(param_vals_baseline, param_names, modParam_scaling, modParam_names, nTrials)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

param_vals_scaled = repmat(param_vals_baseline.', nTrials, 1);
modParam_baseline = zeros(numel(modParam_names), 1);
modParam_vals = zeros(nTrials, numel(modParam_names));

    for mod_index = 1:(length(modParam_names))
        param_name = modParam_names(mod_index);
        index = find(strcmp(param_names, param_name));
        modParam_baseline(mod_index) = param_vals_baseline(index);
        modParam_vals(:, mod_index) = param_vals_baseline(index) * modParam_scaling(:,mod_index);
        param_vals_scaled(:, index) = modParam_vals(:, mod_index);
    end

end