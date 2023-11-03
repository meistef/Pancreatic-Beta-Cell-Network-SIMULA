function B = set_up_matrix(G)
% B = set_up_matrix(G) Set up matrix for the discrete model

% Load parameters
N = G.N;
Cm = G.Cm;
dt = G.dt;
Am = reshape(G.Am, G.N, 1);

% Set up sub matrices
I = speye(N,N);
Bi = set_up_sub_matrix(G, G.Mi);

% Combine matrices
B = I+(dt/Cm)*Bi./Am;

B = sparse(B);
end
