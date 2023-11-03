function G = read_islet_data(basename, exp_name, only_beta)
%G = read_islet_data(basename, exp_name)

% Read geometric information
filename = sprintf('%s/%s/%s_Reconstructed_Arch.txt', basename, exp_name, exp_name);
fileID = fopen(filename, 'r');
formatSpec = '%f %f %f %f %f %f';
sizeA = [6, inf];
A = fscanf(fileID, formatSpec, sizeA);
fclose(fileID);
A = A';
G.R = A(:,1);
G.types = round(A(:,3));
G.types = G.types - 11; % 0 is alpha, 1 is beta
G.coords = A(:,4:6);
G.N = length(G.R);

% Read network information
filename = sprintf('%s/%s/%s_global_contacts.txt', basename, exp_name, exp_name);
fileID = fopen(filename, 'r');
formatSpec = '%f';
B = fscanf(fileID, formatSpec);
fclose(fileID);
B = reshape(B, size(A, 1), size(A, 1));
G.contacts = B;

if only_beta
    include_idx = find(G.types == 1);
    G.R = G.R(include_idx);
    G.types = G.types(include_idx);
    G.coords = G.coords(include_idx, :);
    G.N = length(include_idx);
    G.contacts = G.contacts(include_idx,:);
    G.contacts = G.contacts(:, include_idx);
end

end