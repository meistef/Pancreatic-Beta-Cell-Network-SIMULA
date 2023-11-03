function [V_m, states_m] = solve_system(G)
% [V_m, states_m] = solve_system(G)
% Run a simulation of the monodomain model specified by the parameters G


% Set up matrix
A = set_up_matrix(G);
if strcmp(G.solution_method, 'iterative')
    [LA, UA] = ilu(A);
end

% Set up default initial conditions and parameters
states = G.init_states();
Ns = length(states);
param = G.init_parameters();

% Set up intial conditions for the entire domain
states = states*ones(1, G.N);
t = 0;
V = states(G.V_idx,:)';

% Set up cell model parameters
P = set_up_param_from_file(param, G);

% Set up zero vector
d = zeros(G.N, 1);

if isfield(G, 'DT')
    G.num_save = round(G.Tstop/G.DT) + 1;
    G.save_step = round(G.DT/G.dt);
else
    G.num_save = G.Nt+1;
    G.save_step = 1;
    G.DT = G.dt;
end


% Matrix for saving the solution
V_m = zeros(G.N, G.num_save);
V_m(:,1) = V;
states_m = zeros(G.N, Ns, G.num_save);
states_m(:,:,1) = states';

% Run simulation
fprintf('Starting simulation...\n')
print_step = 500*G.dt;
n_print = max(round(print_step/G.dt), 1);
t1 = tic;
for n = 1:G.Nt
    
    % Step 1: Solve ode-system for membrane model (forward Euler)
    for k=1:G.nt
        states = states + G.dt_ode*G.rhs(t, states, P);
        t = t + G.dt_ode;
    end
    V = states(G.V_idx, :)';
    if any(isnan(V)) || any(isinf(V))
        fprintf('ODE solution is inf or nan.\n')
        break;
    end
    
    
    % Step 2: Solve pde-system
    b = V;
    x0 = V;
    if strcmp(G.solution_method, 'iterative')
        [X, ~] = bicgstab(A, b, G.tol, 1000, LA, UA, x0);
    else
        X = A\b;
    end
    
    % Extract and save the solution
    V = X(1:G.N);
    
    if rem(n, G.save_step) == 0
        V_m(:, round(n/G.save_step)+1) = V;
        states_m(:,:,round(n/G.save_step)+1) = states';
    end
   
    % Update the membrane potential in the cell model
    states(G.V_idx, :) = V';
    
    % Estimate remaining simulation time
    if rem(n, n_print) == 0
        % Print current point in time
        fprintf('t = %.2f sec. ', t/1000);

        % Print estimated simulation time
        t2 = toc(t1);                % Time usage for n_print time steps
        t_rem = t2*(G.Nt-n)/n_print; % Estimated remaining simulation time
        fprintf('Estimated remaining simulation time: ');
        if t_rem > 86400*2
            fprintf('%.1f days \n', t_rem/86400);
        elseif t_rem > 3600
            fprintf('%.1f h \n', t_rem/3600);
        elseif t_rem > 60
            fprintf('%.1f min \n', t_rem/60);
        else
            fprintf('%.1f sec \n', t_rem);
        end
        t1 = tic;
        
    end
end
fprintf('Simulation done\n\n')



end


