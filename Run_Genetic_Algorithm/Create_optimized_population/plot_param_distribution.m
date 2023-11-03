% Code last modified: Roshni (11 August 2023)
% Comments added 2 November 2023

clear; clc; close all

currentDir = pwd;
addpath(genpath(currentDir(1:find(currentDir==filesep,1,'last')-1))); % add all subfolders of the working folder's parent folder (this adds all model-specific folders to Matlab's available directories)

%%
%Final population stdevs
load('modParam_scaling_1.00.mat')
modParam_names = {'V_{GKmax}', 'K_{GK}', 'V_{PFKmax}', 'K_{PFK}', 'h_{PFK}', 'V_{GAPDHmax}', 'g_{Kv}', 'g_{BK}', 'g_{CaL}', 'g_{CaPQ}', 'g_{CaT}', 'g_{KATP}', 'g_{HERG}', 'V_{hNa}', 'n_{hNa}', 'g_{Na}', 'V_{mNa}', 'V_{hNalow}','n_{hNalow}', 'g_{Nalow}', 'V_{mNalow}'};
modParam_names = convertCharsToStrings(modParam_names);
nModParams = 21;

figure(1); set(gcf, 'color', 'w'); hold on;
for iModParam = 1:nModParams
        subplot(4, ceil(nModParams/4), iModParam);
        iScalingDistribution = modParam_scaling(:,iModParam);
        plotDist = iScalingDistribution>0;
        histogram(iScalingDistribution(plotDist), 'FaceColor', 'blue');
        title(sprintf(' %s: \x03BC = %.2f; \x03C3 = %.2f', convertCharsToStrings(modParam_names{1,iModParam}), mean(iScalingDistribution(plotDist)), std(iScalingDistribution(plotDist))));
        grid on;
end
sgtitle('Scaling Factor Distributions', 'Interpreter', 'none');
figure(1); set(gcf, 'Units', 'Inches', 'Position', [0, 0, 18, 18], 'PaperUnits', 'Inches', 'PaperSize', [18, 18])
%f = gcf; exportgraphics(f,'distribution.png','Resolution', 300)
