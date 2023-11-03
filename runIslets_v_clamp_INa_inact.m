close all; clear all; clc; 
currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

%% Simulation output and simulation parameters
% Fraction of cells expressing high voltage INa
frac_high = 0.15;   %15 percent INa high
frac_low = 1-frac_high; %85% INa low

%saves the output of the ode solver in this uniquely named folder
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1)));
output_general = [pwd,'/Experiments',filesep,'INa_inact_prelim',filesep,'Output',filesep,datestr(now,'yyyymmddTHHMMSS'),filesep,num2str(frac_high),'_',num2str(frac_low)];
if exist(output_general, 'dir') ~= 7
    mkdir(output_general);
end

%% Quick-Access Analysis Settings

%% Baseline PoM parameter settings (no inclusion of INa_inact parameters)
%Check value of G; level of glucose in Riz2014_init_parameters.
nTrials = 3000;          % number of different parameter combinations 
%modParam_stdev = [0.71 0.90 0.03 0.09 0.14 0.27]; % Glycolytic
% modParam_stdev = [0.37 0.80 0.49 0.34 0.40 0.42 0.89 0.30]; % Ionic
modParam_stdev = [0.71 0.90 0.03 0.09 0.14 0.27 0.37 0.80 0.49 0.34 0.40 0.42 0.89 0.30]; % Combined
%modParam_names = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK','V_GAPDH_max'}; % Glycolytic
% modParam_names = {'g_Kv', 'g_BK', 'g_Na', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG'}; % Ionic
modParam_names = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK', 'V_GAPDH_max', 'g_Kv', 'g_BK', 'g_Na', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG'}; % Combined
nModParams = length(modParam_names);
              
%% Define Models to Analyze
% model = [model, {{parameter_function_handle, state_function_handle, rhs_function_handle, 'modelName'}}]; % example format (comment-out when running code)

% Individual models (code lines) can be commented-out without affecting how the code runs/works
model = cell(0);
model = [model, {{@Riz2014_init_parameters_INa_low, @Riz2014_init_states_INa_low, @Riz2014_rhs_INa_low, 'Beta cell'}}];
nModels = numel(model);

%% Setup Model Parameters to be Perturbed

scalingFactorFile = strcat(pwd,'/PoM/PoM_Scaling_Factors/','SAVED_scalingFactors_original.mat'); % Baseline Scaling Factors file not including INa_inact parameter variation
scalingFactorFile_INainact = strcat(pwd,'/PoM/PoM_Scaling_Factors/','SAVED_scalingFactors_INainact.mat'); % Scaling Factors file including INa_inact parameter variation
if isfile(scalingFactorFile_INainact)
    load(scalingFactorFile_INainact);
    fprintf("Full scaling factors loaded.\n");
elseif isfile(scalingFactorFile)
    load(scalingFactorFile);
    fprintf("Baseline scaling factors loaded.\n");
    modParam_scaling(:,find(strcmp(modParam_names,'g_Na'))) = [];
    modParam_stdev(:,find(strcmp(modParam_names,'g_Na'))) = [];
    modParam_names(find(strcmp(modParam_names,'g_Na'))) = [];
    nModParams = nModParams-1;
    INa_high_modParam_names = {'V_hNa', 'n_hNa', 'g_Na', 'V_mNa'};
    INa_high_modParam_stdev = [0.15, 0.15, 0.49, 0.15];
    INa_high_nModParams = numel(INa_high_modParam_names);
    INa_high_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_high_modParam_names, frac_high*nTrials);
    INa_high_modParam_scaling = [INa_high_modParam_scaling;zeros(frac_low*nTrials,INa_high_nModParams)];
    INa_low_modParam_names = {'V_hNa_low','n_hNa_low', 'g_Na_low', 'V_mNa_low'};
    INa_low_modParam_stdev = [0.15, 0.15, 0.49, 0.15];
    INa_low_nModParams = numel(INa_low_modParam_names);
    INa_low_modParam_scaling = getScalingFactors_INa_inact(INa_high_modParam_stdev, INa_low_modParam_names, frac_low*nTrials);
    INa_low_modParam_scaling = [zeros(frac_high*nTrials,INa_low_nModParams);INa_low_modParam_scaling];
    modParam_scaling = [modParam_scaling,INa_high_modParam_scaling,INa_low_modParam_scaling];
    modParam_names = [modParam_names,INa_high_modParam_names, INa_low_modParam_names];
    modParam_stdev = [modParam_stdev,INa_high_modParam_stdev,INa_low_modParam_stdev];
    nModParams = nModParams+INa_high_nModParams+INa_low_nModParams;
    save(scalingFactorFile_INainact, 'modParam_scaling', 'modParam_names', 'nModParams', 'modParam_stdev', 'nTrials', 'frac_high','-v7.3');
end

%%%%%%%%%%%%%%%%

figureFolder_general = [output_general, filesep, 'Figures'];
if exist(figureFolder_general, 'dir') ~= 7
    mkdir(figureFolder_general);
end
plotScalingFactors(modParam_scaling, modParam_names, nModParams, figureFolder_general);

%% 

modParam_vals = cell(nModels, 1);
IC_vals = cell(nModels, 1);
IC_names = cell(nModels, 1);
trial_flags = cell(nModels, 1);

param_vals_scaled = cell(nModels,1); 
modParam_baseline = cell(nModels,1); 

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
 %state_vals_baseline has all the inital conditions for 15 islet
 %params
 nStates = numel(state_vals_baseline);

 
 % if isempty(gcp('nocreate')) % No parallel pool running 
 %     parpool; % Startup parallel processing pool
 % end
   
 % Iteratively solve the model for each parameter trial and store ICs
 runTime = tic;
 fprintf("Running simulation\n");
 %parfor iTrial = 1:10:nTrials 
 parfor iTrial = 1:nTrials  
%  for iTrial = 1:10 %nTrials 
     state_vals_output = state_vals_baseline';
     TrialParamVals = param_vals_scaled(iTrial,:);
     [EE, LE, TE, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa, ramp_IK, Peak_IK_minus20, Peak_IK_max, v_half_act_IK, n_act_IK] = run_v_clamp(param_vals_scaled(iTrial,:), param_names,...
         state_vals_output, state_names, model_rhsFun);
%      regular for loop
%       saveoutputfilepath = [output_general, filesep, num2str(iTrial)];
%       save(saveoutputfilepath, 'EE','LE','TE','v_half', 'vh_r2', 'peak_INa','v_half_act', 'vha_r2', 'peak_ICa','late_ICa','ramp_IK', 'Peak_IK_minus20', 'Peak_IK_max', 'v_half_act_IK', 'n_act_IK', 'TrialParamVals','param_names');
%      parfor loop
     parsave(output_general, EE, LE, TE, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa, ramp_IK, Peak_IK_minus20, Peak_IK_max, v_half_act_IK, n_act_IK, TrialParamVals, param_names, iTrial);
 end


fprintf("\n%.3fs\n", toc(runTime));
clear all;

function parsave(output_general, EE, LE, TE, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa, ramp_IK, Peak_IK_minus20, Peak_IK_max, v_half_act_IK, n_act_IK, TrialParamVals, param_names, i)
  saveoutputfilepath = [output_general, filesep, num2str(i)];
  save(saveoutputfilepath,'EE','LE','TE','v_half','vh_r2', 'peak_INa', 'v_half_act', 'vha_r2', 'peak_ICa','late_ICa','ramp_IK', 'Peak_IK_minus20', 'Peak_IK_max', 'v_half_act_IK', 'n_act_IK','TrialParamVals','param_names');
end






