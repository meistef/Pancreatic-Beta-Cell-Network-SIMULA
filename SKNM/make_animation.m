function make_animation(V, G, filename, unit)
% make_animation(V, G, filename, unit)
% Make an animation of a variable from a discrete model simulation
% 
% Input arguments:
%   V: Solution vector containing the variable to plot
%   G: Domain geometry parameter
%   filename: filename for saving the movie
%   unit: unit of the variable to plot

% Time steps
if ~isfield(G, 'DT')
    G.DT = G.dt;
end
DT = max(10, G.DT);
t = 0:DT:G.Tstop;
t_idx = round(t/G.DT) + 1;

% Maximum and minimum values
V_min = min(min(V));
V_max = max(max(V));

% Set up figure
figure('Units','centimeters', 'Position', [10 10 18 15], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [18, 15]);
set(gcf, 'Color','white')
set(gca, 'FontSize', 18, 'nextplot', 'replacechildren')

% Set up movie writer
writerObj = VideoWriter(filename, 'MPEG-4');
writerObj.FrameRate = 10;
open(writerObj)

% Set up mesh of points
x = G.coords(:,1);
y = G.coords(:,2);
z = G.coords(:,3);

% Axis limits
x_min = min(x-2*G.R);
x_max = max(x+2*G.R);
y_min = min(y-2*G.R);
y_max = max(y+2*G.R);
z_min = min(z-2*G.R);
z_max = max(z+2*G.R);

% Set up plotting parameters
r = 2*G.R;      % Cell diameter for plotting
c_map = jet;  % Colormap for plotting
nc = length(c_map);

% Plot solution
for n=1:length(t_idx)
    
    if n==1
        % Find correct scaling of scatter circles
        h = scatter3(x, y, z, [], c_map(1 + round((nc-1)*(V(:,t_idx(n))-V_min)/(V_max-V_min)),:), 'filled');
        currentunits = get(gca,'Units');
        set(gca, 'Units', 'Points');
        axpos = get(gca,'Position');
        set(gca, 'Units', currentunits);
        markerWidth = r/diff(xlim)*axpos(3);
        set(h, 'SizeData', markerWidth.^2)
    else
        h = scatter3(x, y, z, markerWidth.^2, c_map(1 + round((nc-1)*(V(:,t_idx(n))-V_min)/(V_max-V_min)),:), 'filled');
    end
    
    xlim([x_min, x_max]) 
    ylim([y_min, y_max]) 
    zlim([z_min, z_max]) 

    colormap(c_map)
    c = colorbar;
    caxis([V_min, V_max])
    ylabel(c, unit)
    grid off
    box off
    axis off
    axis equal
    set(gca, 'fontsize', 14)
    title(sprintf('T = %.1f sec', t(n)/1000), 'fontsize', 24)

    % Save frame
    writeVideo(writerObj, getframe(gcf));
end

% Close the video writer
close(writerObj)
end

