function [t, relres, iter] = newton(c, zr, t0, maxiter, tol)
% Solve p(t)=zr for Legendre approx with coeffs c

% Break if diverging
Rbreak = 100; % This is really far from [-1, 1]

% Setup and run Newton
t = t0;
lmax = numel(c)-1;
relres = tol;
for iter=1:maxiter
    [P, D] = legendre.deriv_scalar(lmax, t);
    f = P*c - zr;
    fp = D*c;
    dt = -f / fp;
    t = t + dt;
    relres = abs(dt/t);
    if relres < tol || abs(t) > Rbreak
        break
    end
end

if abs(t) > Rbreak
    % Broke on divergence check, return ugly-looking results
    iter = maxiter;
    relres = 1;
end
