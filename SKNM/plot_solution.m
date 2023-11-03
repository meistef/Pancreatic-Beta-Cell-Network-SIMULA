 function plot_solution(V, G, cell_idx, name, unit)
%plot_solution(V, G, cell_idx, unit)
%   Plot a variable in a cell in the discrete model as a function of time 
% 
% Input arguments:
%   V: Solution vector containing the variable to plot
%   G: Domain geometry parameter
%   cell_idx: index of the cell 
%   name: name of the variable to plot
%   unit: unit of the variable to plot

% Convert cell idx to global idx
% with_cells = zeros(G.Nm, 1);
% with_cells(cell_idx) = 1;
% all_cells = zeros(G.N, 1);
% all_cells(G.with_cell) = with_cells;
% cell_idx = find(all_cells);

% Set up figure
figure('Units','centimeters', 'Position', [10 10 18 15], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [18, 15]);
T = (0:G.DT:G.Tstop);
plot(T/1000, V(cell_idx,:), 'linewidth', 2)
set(gca, 'fontsize', 12)
xlabel('t (sec)')
ylabel(unit)
title(name)
end