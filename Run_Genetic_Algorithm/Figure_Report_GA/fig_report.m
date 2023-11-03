clear; clc; close all

% Code last modified: Roshni (11 August 2023)
% Comments added 2 November 2023

currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

M = readtable('beta_glucose180.csv');  %Experimental values from the Camunas-Soler (2019) paper https://github.com/jcamunas/patchseq

% Experimental Data
Late_CaL_current = M.LateCa2_Current;
early_ca_current = M.EarlyCa2_Current;
peak_INa_current = M.PeakNa_Current;
half_inact_sodium_current = M.HalfInactivationSodiumCurrent_mV;
total_exo = M.TotalExocitosis;
early_exo = M.EarlyExocytosis;

early_ca_current = early_ca_current(~isnan(early_ca_current));
Late_CaL_current = Late_CaL_current(~isnan(Late_CaL_current));Late_CaL_current = Late_CaL_current(Late_CaL_current<=0);
peak_INa_current = peak_INa_current(~isnan(peak_INa_current));
half_inact_sodium_current = half_inact_sodium_current(~isnan(half_inact_sodium_current));
half_inact_sodium_current = half_inact_sodium_current(half_inact_sodium_current<=0);
total_exo = total_exo(~isnan(total_exo)); total_exo = total_exo(total_exo>=0);
early_exo = early_exo(~isnan(early_exo)); early_exo = early_exo(early_exo>=0);

% Simulated population
peak_INa_sim = load('peak_INa_1.00.mat');
early_CaL_sim = load('peak_ICa_1.00.mat');
late_CaL_sim = load('late_ICa_1.00.mat');
v_half_sim = load('v_half_1.00.mat');
ee_sim = load('EE_1.00.mat');
te_sim = load('TE_1.00.mat');
load('modParam_stdev_1.00.mat');
load('frac_high_1.00.mat')

early_CaL_sim = cell2mat(struct2cell(early_CaL_sim));
late_CaL_sim = cell2mat(struct2cell(late_CaL_sim));
peak_INa_sim = cell2mat(struct2cell(peak_INa_sim));
v_half_sim = cell2mat(struct2cell(v_half_sim));
ee_sim = cell2mat(struct2cell(ee_sim));
te_sim = cell2mat(struct2cell(te_sim));

%%
% Normalize data to mean

figure(1); set(gcf, 'color', 'w'); hold on;

subplot(2, 3, 1); hold on
range_earlyCaL = abs(min(early_CaL_sim./mean(early_CaL_sim)) - max(early_CaL_sim./mean(early_CaL_sim)));
bw1 = range_earlyCaL/19;
h1 = histogram(early_ca_current./mean(early_ca_current),'BinWidth',bw1, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
h2 = histogram(early_CaL_sim./mean(early_CaL_sim),'BinWidth',bw1,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Early {\it{Ca^{2+}}} peak'); 
xlabel('Normalized Early {\it{Ca^{2+}}} peak '); ylabel('Frequency')

subplot(2, 3, 2); hold on
range_lateCaL = abs(min(late_CaL_sim./mean(late_CaL_sim)) - max(late_CaL_sim./mean(late_CaL_sim)));
bw2 = range_lateCaL/19;
histogram(Late_CaL_current./mean(Late_CaL_current),'BinWidth',bw2, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(late_CaL_sim./mean(late_CaL_sim),'BinWidth',bw2,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Late {\it{Ca^{2+}}} peak'); 
xlabel('Normalized Late {\it{Ca^{2+}}} peak '); ylabel('Frequency')


subplot(2, 3, 3); hold on
range_peakINa = abs(min(peak_INa_sim./mean(peak_INa_sim)) - max(peak_INa_sim./mean(peak_INa_sim)));
bw3 = range_peakINa/19;
histogram(peak_INa_current./mean(peak_INa_current),'BinWidth',bw3, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(peak_INa_sim./mean(peak_INa_sim),'BinWidth',bw3,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('{\it{I_{Na}}} peak'); 
xlabel('Normalized {\it{I_{Na}}} peak'); ylabel('Frequency')


subplot(2, 3, 4); hold on
range_vhalf = abs(min(half_inact_sodium_current) - max(half_inact_sodium_current));
bw4 = range_vhalf/19;
histogram(half_inact_sodium_current,'BinWidth',bw4, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(v_half_sim,'BinWidth',bw4,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title(' {\it{I_{Na}}} half inactivation'); 
xlabel('{\it{I_{Na}}} half inactivation (mV)'); ylabel('Frequency')

subplot(2, 3, 5); hold on
range_ee = abs(min(early_exo./mean(early_exo)) - max(early_exo./mean(early_exo)));
bw5 = range_ee/19;
histogram(early_exo./mean(early_exo),'BinWidth',bw5, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(ee_sim./mean(ee_sim),'BinWidth',bw5,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Early Exocytosis'); 
xlabel('Normalized Early Exocytosis'); ylabel('Frequency')

subplot(2, 3, 6); hold on
range_te = abs(min(total_exo./mean(total_exo)) - max(total_exo./mean(total_exo)));
bw6 = range_te/19;
histogram(total_exo./mean(total_exo),'BinWidth',bw6, 'Normalization', 'Probability','facecolor','b', 'facealpha',.4);
hold on; 
histogram(te_sim./mean(te_sim),'BinWidth',bw6,'Normalization', 'Probability', 'facecolor','r','facealpha',.4);
title('Total Exocytosis'); 
xlabel('Normalized Total Exocytosis'); ylabel('Frequency')
legend('Experimental Distribution', 'Simulated Distribution','NumColumns',2)

%%
set(findobj(gcf,'type','axes'), 'FontName','Arial','FontSize',13.5, 'LineWidth', 1, 'box', 'off', 'tickdir', 'out');
figure(1); set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10, 8])
%f = gcf; exportgraphics(f,'report_fig_Roshni.png','Resolution', 300)


