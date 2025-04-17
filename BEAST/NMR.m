% --------------------------------------------------------------
% This script defines the parameters for a dynamical system that 
% models NMR (Nuclear Magnetic Resonance) applications.
% --------------------------------------------------------------

close all

% ------------------------------------------------------------------
% Initialization Flags
% ------------------------------------------------------------------
initialize = true;         % Initializes system with deterministic solve
Schrodinger = false;       % Vectorizes state into real and imaginary parts
min_energy = true;         % Minimizes energy
plot_sol = true;           % Plots the solution
report_progress = true;    % Reports progress

% ------------------------------------------------------------------
% Time Parameters
% ------------------------------------------------------------------
par.T = pi;                % Total time duration
par.K = 150;               % Number of time steps

% ------------------------------------------------------------------
% Moment Polynomial Order and Uncertainty Spreads
% ------------------------------------------------------------------
par.Na = 1;                % Moment polynomial order for system A
par.Nb = 1;                % Moment polynomial order for system B
par.alpha_min = -1;        % Minimum uncertainty for alpha
par.alpha_max = 1;         % Maximum uncertainty for alpha
par.beta_min = 0.9;        % Minimum uncertainty for beta
par.beta_max = 1.1;        % Maximum uncertainty for beta

% ------------------------------------------------------------------
% Adjustable System Parameters
% ------------------------------------------------------------------
par.epsilon = 3e-3;        % Epsilon for convergence criteria
par.iter_max = 150;        % Maximum number of iterations
par.lambda = 0.01;         % Regularization parameter (lambda)
par.mu = 7 * par.Na * par.Nb;  % Regularization parameter (mu)

% ------------------------------------------------------------------
% Initial and Target States
% ------------------------------------------------------------------
par.n_base = 3;            % Number of states (real)
par.m = 2;                 % Number of control inputs
par.X0 = [0; 1; 1];        % Initial state vector
par.XT = [1; 1; 0];        % Target state vector
par.n_target = find(par.XT);  % Index of target state

% ------------------------------------------------------------------
% Control Bounds (Input Constraints)
% ------------------------------------------------------------------
par.u0 = zeros(par.m, par.K - 1);  % Initial control inputs (zeros)
par.u_max = Inf;                   % Maximum control input (infinity)
par.u_min = -Inf;                  % Minimum control input (negative infinity)
par.du_max = Inf;                  % Maximum change in control input (infinity)
par.du_min = -Inf;                 % Minimum change in control input (negative infinity)

% ------------------------------------------------------------------
% State Matrices (System Dynamics)
% ------------------------------------------------------------------
% Base matrices for A, B1, and B2 that describe the system dynamics
par.A_base = spalloc(3, 3, 2);      % Sparse matrix for system A
par.w = 1;                          % Weight factor for A matrix
par.A_base(2, 1) = -par.w;          % Set A matrix values
par.A_base(1, 2) = par.w;           % Set A matrix values

% Define control matrices B1 and B2
par.B1_base = spalloc(3, 3, 2);     % Sparse matrix for B1
par.B1_base(1, 3) = 1;              % Set B1 matrix values
par.B1_base(3, 1) = -1;             % Set B1 matrix values

par.B2_base = spalloc(3, 3, 2);     % Sparse matrix for B2
par.B2_base(2, 3) = -1;             % Set B2 matrix values
par.B2_base(3, 2) = 1;              % Set B2 matrix values

% Combine B1 and B2 into a single matrix B_base
par.B_base = [par.B1_base, par.B2_base];  % Concatenate B1 and B2 matrices

% ------------------------------------------------------------------
% Start Timer for Performance Measurement
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