% Run an example simulation of SKNM for beta cells
clear all
restoredefaultpath;
addpath('Riz_2014_INa_low/') % Membrane model

% Set up simulation time
Tstop = 2e3; % Total simulation time (ms)

% Read input data
only_beta = 1;
basename = 'Human/';
exp_name = 'H06';
extra_name = '_example';
filename = sprintf('%s%s', exp_name, extra_name); % Filename for animation
G = read_islet_data(basename, exp_name, only_beta);

% Shift indices
[G, G_old] = shift_indices(G);

% Generate model parameters
G = domain_geometry(G, Tstop);
G.Mi = set_up_conductivity(G);

% Specify which parameters should be varied between cells and specify files 
% in which the parameter values for the individual cells are specified
G.parameters_from_file = {'V_GK_max', 'K_GK', 'V_PFK_max', 'K_PFK', 'h_PFK', ...
    'V_GAPDH_max', 'g_Kv', 'g_BK', 'g_CaL', 'g_CaPQ', 'g_CaT', 'g_KATP_hat', 'g_HERG', ...
    'V_hNa', 'n_hNa', 'g_Na', 'V_mNa', 'V_hNa_low', 'n_hNa_low', 'g_Na_low', 'V_mNa_low'};
% find files
% folder_path = "cell_property_distributions/hg/";
% filelist = dir(fullfile(folder_path, '*.txt*'));
% filenames = {};
% for i = 1:length(filelist)
%     filenames{end+1} = filelist(i).name;
% end

G.parameter_files = {'cell_property_distributions/hg/V_GK_max_hg.txt', ...
    'cell_property_distributions/hg/K_GK_hg.txt', ...
    'cell_property_distributions/hg/V_PFK_max_hg.txt', ...
    'cell_property_distributions/hg/K_PFK_hg.txt', ...
    'cell_property_distributions/hg/h_PFK_hg.txt', ...
    'cell_property_distributions/hg/V_GAPDH_max_hg.txt', ...
    'cell_property_distributions/hg/g_Kv_hg.txt', ...
    'cell_property_distributions/hg/g_BK_hg.txt', ...
    'cell_property_distributions/hg/g_CaL_hg.txt', ...
    'cell_property_distributions/hg/g_CaPQ_hg.txt', ...
    'cell_property_distributions/hg/g_CaT_hg.txt', ...
    'cell_property_distributions/hg/g_KATP_hat_hg.txt', ...
    'cell_property_distributions/hg/g_HERG_hg.txt', ...
    'cell_property_distributions/hg/V_hNa_hg.txt', ...
    'cell_property_distributions/hg/n_hNa_hg.txt', ...
    'cell_property_distributions/hg/g_Na_hg.txt', ...
    'cell_property_distributions/hg/V_mNa_hg.txt', ...
    'cell_property_distributions/hg/V_hNa_low_hg.txt', ...
    'cell_property_distributions/hg/n_hNa_low_hg.txt', ...
    'cell_property_distributions/hg/g_Na_low_hg.txt', ...
    'cell_property_distributions/hg/V_mNa_low_hg.txt',};


% Run a simulation
% Here, V contains the membrane potential and states contains all state 
% variables in the membrane model
tic
[V, states] = solve_system(G);
toc

% Extract calcium concentrations
Ca_m = reshape(states(:,G.Ca_m_idx,:), size(V));
Ca_c = reshape(states(:,G.Ca_c_idx,:), size(V));


% Make animations and plots
make_animation(V, G, sprintf('movies/%s_V.mp4', filename), 'mV')
make_animation(Ca_c, G, sprintf('movies/%s_Ca_c.mp4', filename), '\muM')

plot_solution(V, G, 100, 'v', 'mV')
