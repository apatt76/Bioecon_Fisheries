function eta = generate_OU_noise(n_species, n_steps, tau, sigma, dt, rho)
% Generate pre-computed OU noise traces with shared + independent drivers
%
% Inputs:
%   n_species : number of species
%   n_steps   : number of simulation timesteps
%   tau       : autocorrelation timescale (same units as dt)
%   sigma     : total noise amplitude (stationary std dev)
%   dt        : simulation timestep
%   rho       : fraction of variance from shared driver (0 = fully
%               independent, 1 = fully shared). Default 0 if omitted.
%
% Output:
%   eta       : (n_species x n_steps) matrix of noise values

if nargin < 6
    rho = 0;
end

phi        = exp(-dt / tau);
sigma_step = sigma * sqrt(1 - phi^2);   % innovation SD for target stationary variance

% Shared and independent innovation SDs scale by sqrt(rho) and sqrt(1-rho)
% so that variances (not SDs) partition correctly
sigma_shared = sigma_step * sqrt(rho);
sigma_indep  = sigma_step * sqrt(1 - rho);

% Draw all innovations at once — one shared trace, one per-species matrix
innov_shared = sigma_shared * randn(1, n_steps);
innov_indep  = sigma_indep  * randn(n_species, n_steps);

% Initialize shared and independent OU states from stationary distribution
sigma_stat   = sigma;                    % stationary SD of total eta
eta_shared   = sqrt(rho)       * sigma_stat * randn(1, 1);
eta_indep    = sqrt(1 - rho)   * sigma_stat * randn(n_species, 1);

eta = zeros(n_species, n_steps);
eta(:, 1) = eta_shared + eta_indep;

for t = 2:n_steps
    eta_shared = phi * eta_shared + innov_shared(:, t);
    eta_indep  = phi * eta_indep  + innov_indep(:, t);
    eta(:, t)  = eta_shared + eta_indep;
end