%% Solid-Liquid Extraction (SLE) Analytical Simulation

% Computes and plots the volume-average concentrations, extraction yield, 
% and spatiotemporal concentration profile in the solid phase.

clc;
clear;
close all;

%% 1. Model Parameters & Dimensionless Numbers

% Number of series terms in the analytical solution, N
N = 10000;
% Biot number (ratio of internal to external mass transfer resistance), Bi
Bi = 10;
% Liquid-to-solid volume ratio, α (a)
a = 2;
% Effective shape dimension, φ (phi)
phi = 1.5;
% Order of the Bessel functions, v
v = (phi-2)/2;
% Ratio of diffusion time to liquid residence time, Φ (Theta)
Theta = 0;
% Equilibrium partition constant, K
K = 1.5;
% Dimensionless inlet concentration of the liquid phase, Win
Win = 0.3;
% Dimensionless initial concentration of the liquid phase, w0
w0 = 0.5;

%% 2. Space & Time Discretization

% Final dimensionless simulation time, τf (tauf)
tauf = 3;
% Number of temporal grid points, Ntau
Ntau = 1000;
% Dimensionless time vector
tau = [linspace(0, tauf, Ntau)]';

% Number of spatial grid points, Nxi
Nxi = 101;
% Dimensionless spatial coordinate vector
xi = [linspace(0, 1, Nxi)]';

%% 3. Analytical Solution Computation

% Compute eigenvalues, λn (Ln), and succesive differences, λn+1-λn (dLn)
[Ln, dLn] = SLE_eigenvalues(N, Bi, a, v, Theta);

% Compute volume-average solid concentration (up), liquid concentration ...
% (w), and yield (eta)
[up, w, eta] = SLE_average_solution(N, Bi, a, v, Theta, K, ...
    Win, w0, Ln, tau);

% Preallocate matrix and compute local solid concentration u(xi, tau) ...
% across all spatial nodes
u = zeros(Ntau,Nxi);
for j = 1:Nxi
    u(:,j) = SLE_local_solution(N, Bi, a, v, Theta, K, ...
        Win, w0, Ln, xi(j), tau);
end

% Generate 2D coordinate grids for 3D surface plotting
[X,T] = meshgrid(xi,tau);

%% 4. Visualization & Plotting

% Figure 1: Solid Phase Volume-Average Concentration
figure('Name', 'Solid Phase Concentration');
plot(tau, up, 'LineWidth', 1.5);
grid on;
title('Volume-Average Concentration in Solid Phase', ...
    'Interpreter', 'latex');
xlabel('Dimensionless Time, $\tau$', 'Interpreter', 'latex');
ylabel('Average Solid Concentration, $\bar{u}$', 'Interpreter', 'latex');

% Figure 2: Liquid Phase Concentration
figure('Name', 'Liquid Phase Concentration');
plot(tau, w, 'LineWidth', 1.5);
grid on;
title('Volume-Average Concentration in Liquid Phase', ...
    'Interpreter', 'latex');
xlabel('Dimensionless Time, $\tau$', 'Interpreter', 'latex');
ylabel('Liquid Concentration, $w$','Interpreter','latex');

% Figure 3: Extraction Yield
figure('Name', 'Extraction Yield');
plot(tau, eta, 'LineWidth', 1.5);
grid on;
title('Extraction Yield','Interpreter','latex');
xlabel('Dimensionless Time, $\tau$', 'Interpreter', 'latex');
ylabel('Extraction Yield, $\eta$', 'Interpreter', 'latex');

% Figure 4: Spatiotemporal Concentration Profile in Solid Phase
figure('Name', 'Solid Concentration Profile');
surf(X, T, u, 'EdgeColor', 'none');
title('Concentration Profile in Solid Phase', 'Interpreter', 'latex');
xlabel('\xi');
ylabel('Dimensionless Time, $\tau$', 'Interpreter', 'latex');
zlabel('Local Concentration, $u(\xi,\tau)$', 'Interpreter', 'latex');
colorbar;
view(45,30);
axis tight;
axis([0 1 0 1.5 0 1]);