% Code last modified: Roshni (11 August 2023)
% Comments added 2 November 2023

clear; clc; close all

currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

%%

nTrials = 257;  %Number of cells in final population

frac_high = 0.21;  % Optimized using GA
frac_low = round(1-frac_high,2);   

trials_low = round(nTrials*frac_low);
trials_high = nTrials - trials_low;

%5
% modParam_stdev_noNa = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.34 0.40 0.42 0.89 0.30]; 
% INa_high_modParam_stdev = [0.15, 0.15, 0.49, 0.15];
% INa_low_modParam_stdev = [0.15, 0.15, 0.49, 0.15];

%1 Stdev optimized by GA + experimentally know stdev for other parameters
modParam_stdev_noNa = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.156 0.536 0.42 0.89 0.30];
INa_high_modParam_stdev = [0.1762, 0.15, 0.9073, 0.2631];
INa_low_modParam_stdev = [0.2498, 0.15, 0.6788, 0.0547];

% modParam_stdev_noNa = [0 0 0 0 0 0 0 0 0 0.156 0.536 0 0 0]; 
% INa_high_modParam_stdev = [0.1762, 0, 0.9073, 0.2631];
% INa_low_modParam_stdev = [0.2498, 0, 0.6788, 0.0547];

%3
% modParam_stdev_noNa = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.156 0.536 0.42 0.89 0.30]; % standard deviations (normalized to the mean) of the parameter names below
% INa_high_modParam_stdev = [0.1835, 0.15, 0.8528, 0.1479];
% INa_low_modParam_stdev = [0.1912, 0.15, 0.5283, 0.1912];

% modParam_stdev_noNa = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.156 0.536 0.42 0.89 0.30]; % standard deviations (normalized to the mean) of the parameter names below
% INa_high_modParam_stdev = [0.1835, 0.15, 0.4271, 0.1455];
% INa_low_modParam_stdev = [0.0774, 0.15, 0.5283, 0.2850];

%AUG13
% modParam_stdev_noNa = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.156 0.536 0.709 0.89 0.30]; % standard deviations (normalized to the mean) of the parameter names below
% INa_high_modParam_stdev = [0.054, 0.15, 0.4787, 0.0861];
% INa_low_modParam_stdev = [0.2604, 0.15, 0.4376, 0.2604];


modParam_names_noNa = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK', 'V_GAPDH_max', 'g_Kv', 'g_BK', 'g_Na', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG'};
%INa_high_modParam_names = {'V_hNa', 'n_hNa', 'g_Na', 'V_mNa'};
%INa_low_modParam_names = {'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'};
nModParams_noNa = numel(modParam_names_noNa);
modParam_scaling_noNa = getScalingFactors(modParam_stdev_noNa, nModParams_noNa, nTrials);

%% Define Models to Analyze
model = cell(0);
model = [model, {{@Riz2014_init_parameters_INa_low, @Riz2014_init_states_INa_low, @Riz2014_rhs_INa_low, 'Beta cell'}}];


%% INa_high_low
modParam_scaling_noNa(:,find(strcmp(modParam_names_noNa,'g_Na'))) = [];
modParam_stdev_noNa(:,find(strcmp(modParam_names_noNa,'g_Na'))) = [];
modParam_names_noNa(find(strcmp(modParam_names_noNa,'g_Na'))) = [];
nModParams = nModParams_noNa-1;

INa_high_modParam_names = {'V_hNa', 'n_hNa', 'g_Na', 'V_mNa'};
INa_high_nModParams = numel(INa_high_modParam_names);
INa_high_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_high_modParam_names, trials_high);
INa_high_modParam_scaling = [INa_high_modParam_scaling; zeros(trials_low,INa_high_nModParams)];

INa_low_modParam_names = {'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'};
INa_low_nModParams = numel(INa_low_modParam_names);
INa_low_modParam_scaling = getScalingFactors_INa_inact(INa_low_modParam_stdev, INa_low_modParam_names, trials_low);
INa_low_modParam_scaling = [zeros(trials_high,INa_low_nModParams);INa_low_modParam_scaling];

modParam_scaling = [modParam_scaling_noNa,INa_high_modParam_scaling,INa_low_modParam_scaling];
modParam_names = [modParam_names_noNa,INa_high_modParam_names, INa_low_modParam_names];
modParam_stdev = [modParam_stdev_noNa,INa_high_modParam_stdev,INa_low_modParam_stdev];
nModParams = nModParams+INa_high_nModParams+INa_low_nModParams;

% figureFolder_general = [output_general, filesep, 'Figures'];
% if exist(figureFolder_general, 'dir') ~= 7
%     mkdir(figureFolder_general);
% end
% plotScalingFactors(modParam_scaling, modParam_names, nModParams, figureFolder_general);

model_paramFun = model{1}{1};
model_stateFun = model{1}{2};
model_rhsFun   = model{1}{3};
model_name     = model{1}{4};

%%
[param_vals_baseline, param_names] = model_paramFun();
[param_vals_scaled, modParam_baseline, modParam_vals] = modifyParams(...
    param_vals_baseline, param_names, modParam_scaling,...
    modParam_names, nTrials);

% Retrieve Baseline State Values
[state_vals_baseline, state_names] = model_stateFun();
%state_vals_baseline has al l the inital conditions for 15 isletparams
nStates = numel(state_vals_baseline);
EE = zeros(nTrials, 1); TE = zeros(nTrials, 1); 
peak_INa = zeros(nTrials, 1); v_half_act = zeros(nTrials, 1); vha_r2 = zeros(nTrials, 1); peak_ICa = zeros(nTrials, 1); late_ICa = zeros(nTrials, 1);
v_half = zeros(nTrials, 1); vh_r2 = zeros(nTrials, 1);

%Compute output metrics of the population

for iTrial = 1:nTrials 
    state_vals_output = state_vals_baseline';
    
    [EE(iTrial, 1), ~, TE(iTrial,1)] = Exocytosis(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
    [peak_INa(iTrial, 1), v_half_act(iTrial, 1), vha_r2(iTrial, 1), peak_ICa(iTrial, 1), late_ICa(iTrial, 1)] = PeakCurrents(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);
    [v_half(iTrial, 1), vh_r2(iTrial, 1)] = INa_inact(param_vals_scaled(iTrial,:), param_names, state_vals_output, state_names, model_rhsFun);

end

% Save output metrics
num = 1;
scaling_name = sprintf('modParam_scaling_%.2f.mat', num); save(scaling_name, 'modParam_scaling');
stdev_name = sprintf('modParam_stdev_%.2f.mat', num); save(stdev_name, 'modParam_stdev');
peak_INa_name = sprintf('peak_INa_%.2f.mat', num); save(peak_INa_name, 'peak_INa');
v_half_act_name = sprintf('v_half_act_%.2f.mat', num); save(v_half_act_name, 'v_half_act');
peak_ICa_name = sprintf('peak_ICa_%.2f.mat', num); save(peak_ICa_name, 'peak_ICa');
late_ICa_name = sprintf('late_ICa_%.2f.mat', num); save(late_ICa_name, 'late_ICa');
EE_name = sprintf('EE_%.2f.mat', num); save(EE_name, 'EE');
TE_name = sprintf('TE_%.2f.mat', num); save(TE_name, 'TE');
v_half_name = sprintf('v_half_%.2f.mat', num); save(v_half_name, 'v_half');
frac_high_name = sprintf('frac_high_%.2f.mat', num); save(frac_high_name, 'frac_high');

