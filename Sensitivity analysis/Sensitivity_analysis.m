%% Set-up directories and data
% Load the initial parameters
clear all; clc; close all
currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

% Create a cell array of labels for the y-axis
outcome_labels = {'EE', 'LE', 'TE', 'v_half', 'vh_r2', 'peak_INa', 'v_half_act', ...
    'vha_r2', 'peak_ICa', 'late_ICa'};
outcomes_labels_select = {'EE', 'TE', 'v_half', 'peak_INa', 'peak_ICa', 'late_ICa'};

% Initialize to store the data
outcomes_dict_select = cell(1, length(outcome_labels));

% Obtain the baseline parameters and parameter names
[init_params, param_names] = Riz2014_init_parameters_INa_low();
merged_df = [param_names, num2cell(init_params)];

% Load the initial states
[state_vals_baseline, state_names] = Riz2014_init_states_INa_low();
state_vals_output = state_vals_baseline';

% Create a vector that will be used for the sensitivity analysis
sens_scale = [4, 2, 1, 0.5, 0.25];

%% Sensitivity analysis from scratch

% Initialize cell arrays to store the scaled result
EE = zeros(size(merged_df, 1), numel(sens_scale));
LE = EE;
TE = EE;
v_half = EE;
vh_r2 = EE;
peak_INa = EE;
v_half_act = EE;
vha_r2 = EE;
peak_ICa = EE;
late_ICa = EE;

% Loop through the dataframe and scale each parameter 
tic 
for i = 1:size(merged_df, 1)
    fprintf('Parameter #%d', i)
    for j = 1:numel(sens_scale)
        X_scale = cell2mat(merged_df(:, 2)); 
        X_scale(i) = cell2mat(merged_df(i, 2))* sens_scale(j);
        [EE(i, j), LE(i, j), TE(i, j), v_half(i, j), vh_r2(i, j), peak_INa(i, j),...
            v_half_act(i, j), vha_r2(i, j), peak_ICa(i, j), late_ICa(i, j)] = run_v_clamp_sens(X_scale,...
            param_names, state_vals_output, state_names, @Riz2014_rhs_INa_low);
    end
end
toc

% Create a list with all results and selected params
outcomes_dict = {EE, LE, TE, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa};
outcomes_dict_select = {EE, TE, v_half, peak_INa, peak_ICa, late_ICa};

% Create a cell array of labels for the y-axis
outcome_labels = {'EE', 'LE', 'TE', 'v_half', 'vh_r2', 'peak_INa', 'v_half_act', ...
    'vha_r2', 'peak_ICa', 'late_ICa'};
outcome_labels_select = {'EE', 'TE', 'v_half', 'peak_INa', 'peak_ICa', 'late_ICa'};

% Loop through to store outcomes
for i = 1:size(outcomes_dict_select, 2)
    csv_name = sprintf('sens_outcomes_%s.csv', outcome_labels_select{i});
    writematrix(outcomes_dict_select{i}, csv_name);
end

%% Calculate the relative difference and visualize the results

% Initialize an empty matrix to store the columns
four_times = [];
two_times = [];
half_times = [];
quarter_times = [];

% Loop through the outcomes_dict and calculate ratios
ratios_dict = cell(size(outcomes_dict_select));
for i = 1:numel(outcomes_dict_select)
    
    % Subset the matrices
    outcome_matrix = outcomes_dict_select{i};
    
    % Subset the scaling factors
    four_time = outcome_matrix(:, 1);
    two_time = outcome_matrix(:, 2);
    half_time = outcome_matrix(:, 4);
    quarter_time = outcome_matrix(:, 5);

    % Calculate the ratios relative to 1x scaling
    ratios_matrix = outcome_matrix ./ outcome_matrix(:, 3);
    ratios_dict{i}= ratios_matrix;
    
    % Stack the scaled columns into the empty matrices

    four_times = [four_times, round(ratios_dict{i}(:,1),2)];
    two_times = [two_times, round(ratios_dict{i}(:,2),2)];
    half_times = [half_times, round(ratios_dict{i}(:,4),2)];
    quarter_times = [quarter_times, round(ratios_dict{i}(:,5),2)];
    
end

% Combined high and low parameters for indexing purpose
modParam_names = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK', 'V_GAPDH_max', 'g_Kv', 'g_BK', 'g_Na', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG', 'V_hNa', 'n_hNa', 'V_mNa', 'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'}; % Combined

% Loop through to store outcomes
param_idx = [];
for i = 1:size(modParam_names, 2)
    param_idx{i} = find(strcmp(param_names, modParam_names{i}));
end

% Index the matrices
four_times_idx = four_times(cell2mat(param_idx), :);
two_times_idx = two_times(cell2mat(param_idx), :);
half_times_idx = half_times(cell2mat(param_idx), :);
quarter_times_idx = quarter_times(cell2mat(param_idx), :);

% Relabel axes
param_names_y = {'V_{GKmax}', 'K_{GK}', 'V_{PFKmax}', 'K_{PFK}', 'h_{PFK}', 'V_{GAPDHmax}', 'g_{Kv}', 'g_{BK}', 'g_{Na}', 'g_{CaL}', 'g_{CaPQ}', 'g_{CaT}', 'g_{KATP}', 'g_{HERG}', 'V_{hNa}', 'n_{hNa}', 'V_{mNa}', 'V_{hNalow}','n_{hNalow}', 'g_{Nalow}', 'V_{mNalow}'};
outcome_labels_x = {'EE', 'LE', 'TE', 'v_{half}', 'vh_{r2}', 'peak_{INa}', 'v_{halfact}', ...
    'vha_{r2}', 'peak_{ICa}', 'late_{ICa}'};
outcome_labels_x = {'EE', 'TE', 'v_{half}', 'peak_{INa}', 'peak_{ICa}', 'late_{ICa}'};

% EE, TE, vhalf, peakIna, peakIca, lateIca

%% Red to white colormap 

% Create a heatmap for each of scaling factors
figure(1); set(gcf, 'color', 'w')

nColors = 64;
whiteColor = [1 1 1];
redColor = [1 0 0];
blueColor = [0 0 1];

% Linearly interpolate between blue and white (reversed)
blueToWhite = [linspace(blueColor(1), whiteColor(1), nColors/2); ...
               linspace(blueColor(2), whiteColor(2), nColors/2); ...
               linspace(blueColor(3), whiteColor(3), nColors/2)]';

% Linearly interpolate between white and red (reversed)
whiteToRed = [linspace(whiteColor(1), redColor(1), nColors/2); ...
              linspace(whiteColor(2), redColor(2), nColors/2); ...
              linspace(whiteColor(3), redColor(3), nColors/2)]';

redWhiteBlueColormap = [blueToWhite; whiteToRed];
colormap(redWhiteBlueColormap);

subplot(2, 2, 1);
heat_four = heatmap(outcome_labels_x, param_names_y, four_times_idx, 'FontName', 'Arial');
title('4x Scalar');
colormap(redWhiteBlueColormap); caxis([-3.5 5.5]) 

subplot(2, 2, 2);
heat_two = heatmap(outcome_labels_x, param_names_y, two_times_idx, 'FontName', 'Arial');
title('2x Scalar');
colormap(redWhiteBlueColormap); caxis([-1, 3])

subplot(2, 2, 3);
heat_half = heatmap(outcome_labels_x, param_names_y, half_times_idx, 'FontName', 'Arial');
title('0.5x Scalar');
colormap(redWhiteBlueColormap); caxis([0, 2])

subplot(2, 2, 4);
heat_quarter = heatmap(outcome_labels_x, param_names_y, quarter_times_idx, 'FontName', 'Arial');
title('0.25x Scalar');
colormap(redWhiteBlueColormap); caxis([0, 2])

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
exportgraphics(gcf, 'sens_analysis_all4.png', 'Resolution', 300);
%%

% Create a heatmap for two scaling factors
figure(2); set(gcf, 'color', 'w')

a = subplot(2, 1, 1);
set(a, 'Position', [.1, 0.525, .8 , .425]);
heat_four = heatmap(outcome_labels_x, param_names_y, two_times_idx, 'CellLabelColor', 'none', 'FontName', 'Arial', 'FontSize', 12);
title('2-fold scalar');
colormap(redWhiteBlueColormap);caxis([-1, 3])

b = subplot(2, 1, 2);
set(b, 'Position', [.1, .025, .8 , .425])
heat_two = heatmap(outcome_labels_x, param_names_y, half_times_idx, 'CellLabelColor', 'none', 'FontName', 'Arial', 'FontSize', 12);
title('0.5-fold scalar');
colormap(redWhiteBlueColormap); caxis([0, 2])

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 14], 'PaperUnits', 'Inches', 'PaperSize', [10, 14])
exportgraphics(gcf, 'sens_analysis_twoandhalf.png', 'Resolution', 300);

