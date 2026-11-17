function [Ln, dLn] = SLE_eigenvalues(N, Bi, a, v, Theta)
% SLE_EIGENVALUES: Computes eigenvalues of the SLE model
%
% Input arguments:
% N:        Number of series terms in the analytical solution
% Bi:       Biot number (dimensionless)
% a:        Liquid-to-solid volume ratio
% v:        Order of the Bessel functions
% Theta:    Ratio of diffusion time to liquid residence time
%
% Output arguments:
% Ln:       Vector of computed eigenvalues
% dLn:      Differences between consecutive eigenvalues

Ln = zeros(N,1);  % Preallocate eigenvalues

x = [0:0.01:N*4]';
f = fLn(x, Bi, a, v, Theta);
idx = find(f(1:end-1).*f(2:end) < 0);
x_init = [x(idx) x(idx+1)];
x_init = x_init(1:N,:);

for n = 1:N
    Ln(n) = fzero(@(x) fLn(x, Bi, a, v, Theta), x_init(n,:));
end

% Compute spacing between successive eigenvalues
dLn = diff(Ln);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fLn(x, Bi, a, v, Theta)
% FLN: Characteristic equation

Jv  = besselj(v, x);
Jv1 = besselj(v+1, x);
phi = 2*v+2;

f = (Theta - x.^2) .* (Bi * Jv - x .* Jv1) - ...
    phi * Bi * x .* Jv1 / a;

end