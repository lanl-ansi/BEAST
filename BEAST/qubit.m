close all

% --------------------------------------------------------------
% This script defines the parameters of a dynamical system that 
% models local control of the Heisenberg interaction of two 
% interacting qubits.
% Reference: Baker, Luke S., et al. "Robust Quantum Gate Preparation 
% in Open Environments." arXiv preprint arXiv:2410.01161 (2024).
% --------------------------------------------------------------

% ------------------------------------------------------------------
% Initialization Flags
% ------------------------------------------------------------------
initialize = true;         % Initializes system with deterministic solve
Schrodinger = true;        % Vectorizes state into real and imaginary parts
min_energy = true;         % Minimizes energy
plot_sol = true;           % Plots the solution
report_progress = true;    % Reports progress

% ------------------------------------------------------------------
% Time Parameters
% ------------------------------------------------------------------
par.T = pi;                   % Total time duration (pi)
par.K = 150;                  % Number of time steps

% ------------------------------------------------------------------
% Moment Polynomial Order and Uncertainty Spreads
% ------------------------------------------------------------------
par.Na = 2;                   % Moment polynomial order for system A
par.Nb = 2;                   % Moment polynomial order for system B
par.alpha_min = 0.9;          % Minimum uncertainty for alpha
par.alpha_max = 1.1;          % Maximum uncertainty for alpha
par.beta_min = 0.9;           % Minimum uncertainty for beta
par.beta_max = 1.1;           % Maximum uncertainty for beta

% ------------------------------------------------------------------
% Initial and Target States
% ------------------------------------------------------------------
par.n_base = 8;               % Number of states (base)
par.m = 4;                    % Number of control inputs
par.X0 = [1; 0; 0; 0; zeros(4,1)]; % Initial state vector
par.XT = [0; 1; 1; 0; zeros(4,1)]; % Target state vector (Bell state)
par.n_target = find(par.XT);  % Index of target state

% ------------------------------------------------------------------
% Adjustable System Parameters
% ------------------------------------------------------------------
par.epsilon = 4e-3;           % Epsilon for convergence criteria
par.iter_max = 300;           % Maximum number of iterations
par.lambda = 0.1;             % State transfer regularization parameter (lambda)
par.mu = 10;                  % Minimal energy regularization parameter (mu)

% ------------------------------------------------------------------
% Control Bounds (Input Constraints)
% ------------------------------------------------------------------
par.u0 = 0.1 * ones(par.m, par.K - 1); % Initial control inputs
par.u_max = Inf;               % Maximum control input
par.u_min = -Inf;              % Minimum control input
par.du_max = Inf;              % Maximum change in control input
par.du_min = -Inf;             % Minimum change in control input

% ------------------------------------------------------------------
% Paulie Matrices for System Dynamics
% ------------------------------------------------------------------
% Define Pauli matrices for the qubits (2x2 matrices)
sigmax = sparse([0, 1; 1, 0]);   % Pauli X (sigma_x)
sigmay = sparse([0, -1; 1, 0]);  % Pauli Y (sigma_y)
sigmaz = sparse([1, 0; 0, -1]);  % Pauli Z (sigma_z)

% Interaction term J (Heisenberg interaction)
J = 2; 

% Hamiltonian A (Heisenberg interaction between two qubits)
par.A_base = J / 4 * (kron(sigmaz, speye(2)) * kron(speye(2), sigmaz) + ...
                     kron(sigmax, speye(2)) * kron(speye(2), sigmax) + ...
                     kron(sigmay, speye(2)) * kron(speye(2), sigmay));

% Control matrices B1, B2, B3, and B4 (control along Pauli matrices)
par.B1_base = 1 / 2 * kron(sigmax, speye(2));
par.B2_base = 1 / 2 * kron(sigmay, speye(2));
par.B3_base = 1 / 2 * kron(speye(2), sigmax);
par.B4_base = 1 / 2 * kron(speye(2), sigmay);

% ------------------------------------------------------------------
% Kronecker Products and Control Matrix Construction
% ------------------------------------------------------------------
% Dynamics of real and imaginary parts of the wave vector
par.A_base = kron(sparse([0, 1; -1, 0]), par.A_base);   % For interaction terms
par.B1_base = kron(sparse([0, 1; -1, 0]), par.B1_base); % For control input B1
par.B2_base = kron(speye(2), par.B2_base);              % For control input B2
par.B3_base = kron(sparse([0, 1; -1, 0]), par.B3_base); % For control input B3
par.B4_base = kron(speye(2), par.B4_base);              % For control input B4

% Combine control matrices into a single matrix
par.B_base = [par.B1_base, par.B2_base, par.B3_base, par.B4_base];

% ------------------------------------------------------------------
% Start Solver
% ------------------------------------------------------------------
tic                   

% Save the current directory and change to the function directory
main_dir = pwd;                   

cd(['funs', filesep]);            % Change to functions directory

% Call the solver with the defined parameters and flags
[x_bar, u_bar, x, par] = solver(par, initialize, Schrodinger, min_energy, plot_sol, report_progress);

% Return to the original directory after solver finishes
cd(main_dir);

% Measure and display the run time
run_and_plot_time = toc;
fprintf('Execution Time: %.2f seconds\n', run_and_plot_time);