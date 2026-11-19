function [up,w,eta] = SLE_average_solution(N, Bi, a, v, ...
    Theta, K, Win, w0, Ln, tau)
% SLE_AVERAGE_SOLUTION: Computes the volume-averaged solution for SLE
%
% Input arguments:
% N:        Number of series terms in the analytical solution
% Bi:       Biot number (dimensionless)
% a:        Liquid-to-solid volume ratio
% v:        Order of the Bessel functions
% Theta:    Ratio of diffusion time to liquid residence time
% K:        Equilibrium partition constant
% Win:      Dimensionless inlet concentration of the liquid phase
% w0:       Dimensionless initial concentration of the liquid phase
% Ln:       Eigenvalues of the characteristic equation
% tau:      Fourier number (dimensionless time) 

% Output arguments:
% up:       Volume-averaged concentration of the solid phase
% w:        Concentration of the liquid phase
% eta:      Extraction yield

nt  = length(tau);         % Number of time points
up  = zeros(nt,1);         % Preallocate solution vector
w   = zeros(nt,1);         % Preallocate solution vector
eta = zeros(nt,1);         % Preallocate solution vector

% Compute eigenvalues only if not provided
if isempty(Ln)
    Ln = SLE_eigenvalues(N, Bi, a, v, Theta);
end

Jv  = besselj(v,Ln);
Jv1 = besselj(v+1,Ln);
phi = 2*v+2;

EpL = -2*Ln .* (Bi*Jv - Ln .* Jv1) + ...
    (Theta - Ln.^2) .* ((v * Bi ./ Ln - Ln) .* Jv + (v - Bi) .* Jv1) - ...
    (phi * Bi / a) * (Ln .* Jv - v * Jv1);
An = 2*((Theta * Win - w0 * Ln.^2) .* (Bi * Jv - Ln .* Jv1) - ...
    (phi * K * Bi / a) * Ln .* Jv1) ./ (Ln .* EpL);
Bn = Bi * An .* Jv1 ./ (K * (Bi * Jv - Ln .* Jv1) .* Ln);

if Theta == 0
    A0 = (a * w0 + K)/(a+1);
else
    A0 = Win;
end

% Evaluate the series solution
for j = 1:nt
    Tn = exp(-Ln.^2 * tau(j));
    w(j) = sum(An .* Tn) + A0;
    up(j) = phi * sum(Bn .* Tn) + A0 / K;
    eta(j) = (a * Theta / K) * sum(An .* (1 - Tn) ./ Ln.^2);
end

if Theta == 0
    eta = 1 - up;
end

end