function G = domain_geometry(G, Tstop)
% G = domain_geometry(G, Tstop)

% Set up temporal discretization parameters
G.dt = 1;                    % Time step for the pde step 
G.dt_ode = min(0.1, G.dt);   % Time step for the ode step 
G.Tstop = Tstop;             % Total simulation time
G.nt = round(G.dt/G.dt_ode); % Number of ode steps per pde step 
G.Nt = round(G.Tstop/G.dt);  % Number of pde steps 
G.DT = max(1, G.dt);         % Time step for saving the solution

% Define membrane model
G.init_states = @Riz2014_init_states_INa_low;
G.init_parameters = @Riz2014_init_parameters_INa_low;
[~, G.param_names] = G.init_parameters();
[~, state_names] = G.init_states();
G.V_idx = find(strcmp(state_names, 'V'));
G.Ca_m_idx = find(strcmp(state_names, 'Ca_m'));
G.Ca_c_idx = find(strcmp(state_names, 'Ca_c'));
G.rhs = @Riz2014_rhs_INa_low_vectorized;

% Gap junction resistance
G.GJ_areas = set_up_GJ_areas(G);
mean_area = mean(G.GJ_areas(find(G.GJ_areas>0)));
mean_Gg = 2e-7; % mS
G0 = mean_Gg/mean_area;
G.Gg = G0*G.GJ_areas;

% Geometry: Assumed spheres
G.Am = set_up_cell_area(G);
G.Cm = 1;       % uF/cm^2

% Select solver type
% G.solution_method = 'direct';
G.solution_method = 'iterative';
G.tol = 1e-5;

end
