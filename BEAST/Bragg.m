close all

% --------------------------------------------------------------
% This script defines the parameters for the dynamical system model 
% that describes the matter-wave beamsplitting process in atom 
% interferometry. 
% Reference: Baker, Luke S., et al. "Convergence of iterative 
% quadratic programming for robust fixed-endpoint transfer of 
% bilinear systems." arXiv preprint arXiv:2403.18131 (2024).
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
par.T = 2 * pi;            % Total time duration
par.K = 150;               % Number of time steps

% ------------------------------------------------------------------
% Number of moment polynomials and uncertainty intervals
% ------------------------------------------------------------------
par.Na = 3;                % Number of moment polynomials for system A
par.Nb = 2;                % Number of moment polynomials for system B
par.alpha_min = 0.95;      % Minimum uncertainty for alpha
par.alpha_max = 1.05;      % Maximum uncertainty for alpha
par.beta_min = 0.95;       % Minimum uncertainty for beta
par.beta_max = 1.05;       % Maximum uncertainty for beta

% ------------------------------------------------------------------
% Initial and Target States
% ------------------------------------------------------------------
par.n_base = 6;            % Number of states (real + imaginary)
par.m = 1;                 % Number of control inputs
par.X0 = sparse(zeros(par.n_base,1));  % Initial state vector
par.X0(1) = 1;                         % Set the initial state
par.XT = sparse(zeros(par.n_base,1));  % Target state vector
par.XT(2) = 1;                         % Set the target state
par.n_target = find(par.XT);           % Index of target state

% ------------------------------------------------------------------
% Adjustable System Parameters
% ------------------------------------------------------------------
par.epsilon = 3e-3;          % Epsilon for convergence criteria
par.iter_max = 200;          % Maximum number of iterations
par.lambda = 0.01;           % Regularization parameter (state transfer)
par.mu = 10;                 % Regularization parameter (minimal energy)

% ------------------------------------------------------------------
% Control Bounds (Input Constraints)
% ------------------------------------------------------------------
par.u0 = ones(par.m, par.K - 1);    % Initial control inputs
par.u_max = 30;                     % Maximum control input
par.u_min = 0;                      % Minimum control input
par.du_max = Inf;                   % Maximum change in control input
par.du_min = -Inf;                  % Minimum change in control input

% ------------------------------------------------------------------
% State Matrices (System Dynamics)
% ------------------------------------------------------------------
% Base matrices for A and B
par.A_base = sparse(diag(0:2:par.n_base-2).^2); 
par.B_base = 1/2 * (sparse(diag([sqrt(2), ones(1, par.n_base/2-1-1)], 1) + ...
                          diag([sqrt(2), ones(1, par.n_base/2-1-1)], -1)));

% Taking real and imaginary parts
if Schrodinger
par.A_base = kron(sparse([0 1; -1 0]), par.A_base); 
par.B_base = kron(sparse([0 1; -1 0]), par.B_base);
end

% ------------------------------------------------------------------
% Start Timer for Performance Measurement
% ------------------------------------------------------------------
tic                   

% Save the current directory, change to function directory for solver
main_dir = pwd;                   

cd(['funs', filesep]);            % Change to functions directory

% Call the solver with the defined parameters and flags
[x_bar, u_bar, x, par] = solver(par, initialize, Schrodinger, min_energy, plot_sol, report_progress);

% Return to the original directory after solver finishes
cd(main_dir);

% Measure and display the run time
run_time = toc;
fprintf('Execution Time: %.2f seconds\n', run_time);