function [ux] = SLE_local_solution(N, Bi, a, v, Theta, K, ...
    Win, w0, Ln, xi, tau)
% SLE_AVERAGE_SOLUTION: Computes the volume-averaged solution for SLE
%
% Input arguments:
% N:        Number of series terms in the analytical solution
% B:        Biot number (dimensionless)
% a:        Liquid-to-solid volume ratio
% v:        Order of the Bessel functions
% Theta:    Ratio of diffusion time to liquid residence time
% K:        Equilibrium partition constant
% Win:      Dimensionless inlet concentration of the liquid phase
% w0:       Dimensionless initial concentration of the liquid phase
% Ln0: Initial guess for the first eigenvalue
% Fo: Fourier number (dimensionless time)
% Ln: Eigenvalues of the characteristic equation
%
% Output arguments:
% Y: Volume-averaged solution evaluated at Fo

nt = length(tau);          % Number of time points
up = zeros(nt,1);         % Preallocate solution vector
w  = zeros(nt,1);         % Preallocate solution vector

% Compute eigenvalues only if not provided

if isempty(Ln)
    Ln = SLE_eigenvalues(N, Bi, a, v, Theta);
end

Jv  = besselj(v,Ln);
Jvx = besselj(v,Ln*xi);
Jv1 = besselj(v+1,Ln);
phi = 2*v+2;

EpL = -2*Ln .* (Bi*Jv - Ln .* Jv1) + ...
    (Theta - Ln.^2) .* ((v*Bi./Ln - Ln) .* Jv + (v - Bi) .* Jv1) - ...
    (phi*Bi/a) * (Ln .* Jv - v*Jv1);
An = 2*((Theta*Win - w0 * Ln.^2) .* (Bi*Jv - Ln .* Jv1) - ...
    (phi*K*Bi/a)*Ln .* Jv1) ./ (Ln .* EpL);
Bn = Bi * xi.^(-v) * Jvx .* An ./ (K * (Bi*Jv - Ln .* Jv1) );

if Theta == 0
    A0 = (a * w0 + K)/(a+1);
else
    A0 = Win;
end

% Evaluate the series solution
for j = 1:nt
    Tn = exp(-Ln.^2 * tau(j));
    ux(j) = sum(Bn .* Tn) + A0/K;
end

end