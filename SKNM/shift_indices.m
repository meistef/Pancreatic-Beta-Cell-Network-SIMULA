function [G, G_old] = shift_indices(G)

% Input Arguments:
%   G: Domain geometry parameter
%
% Output:
%   G: Domain geometry parameters with updated indices
%   G_old: Initial struct with domain geomoetry parameters
%
% Shift indices of the cells and sort them by their position for
% an easier postprocessing / analysis


x = G.coords(:,1);
y = G.coords(:,2);
z = G.coords(:,3);

% initialize new coordinates and index shift
coords_new = [];
idx_shift = [];
% use mean radius as step size for voxels to search in
step = mean(G.R) * 4;
for zi = min(z):step:max(z)+step
    for yi = min(y):step:max(y)+step
        % get indices for points that are in y and z range
        in_yz = find(yi <= G.coords(:,2) & G.coords(:,2) < yi + step & zi <= G.coords(:,3) & G.coords(:,3) < zi + step);
        % check if points are in voxel
        if isempty(in_yz) == 0
            % get coordinates and append indices
            temp_coords = G.coords(in_yz,:);
            temp_coords = [temp_coords in_yz];
            % sort points with ascending x values
            temp_coords = sortrows(temp_coords, 1);
            % get indices of sorted coordinates
            temp_idx = temp_coords(:,end);

            % append new coordinates and index shifts
            coords_new = [coords_new; temp_coords(:,1:3)];
            idx_shift = [idx_shift; temp_idx];
        end

        in_yz = [];
    end
end

% shift properties to new indices
R_new = G.R(idx_shift);
types_new = G.types(idx_shift);

% initialize new contacts
contacts_new = zeros(length(R_new), length(R_new));
% shift contact values to new indices
for i = 1:length(idx_shift)
    % get row index before shifting
    idx_row_before = idx_shift(i);
    % get nonzero indices in that row
    nonzero = find(G.contacts(idx_row_before, :) ~= 0);
    for j = 1:length(nonzero)
        % get column index of nonzero element
        idx_col_before = nonzero(j);
        % get new column index
        idx_col_now = find(idx_shift == idx_col_before);
        % shift values to new column indices
        contacts_new(i, idx_col_now) = G.contacts(idx_row_before, idx_col_before);
    end
end

% create return values
G_old = G;
G.R = R_new;
G.types = types_new;
G.coords = coords_new;
G.N = G.N;
G.contacts = contacts_new;

end

