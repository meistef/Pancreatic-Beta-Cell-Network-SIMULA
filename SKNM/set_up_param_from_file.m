function P = set_up_param_from_file(param, G)
% P = set_up_param_from_file(param, G)
% Set up parameters for the discrete model with some parameter variations
% specified in files

P = param*ones(1, G.N);

if isfield(G, 'parameters_from_file')
    for n=1:length(G.parameters_from_file)
        param_idx = find(ismember(G.param_names, G.parameters_from_file{n}));
        param_values = importdata(G.parameter_files{n});
        P(param_idx,:) = param_values(1:G.N);
    end
end

end

