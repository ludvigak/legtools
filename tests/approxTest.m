% - function approximation and differentiation
% - quadrature
% - root finding
addpath ..

k = 10;
f = @(t) sin(k*t+1);
fp = @(t) k*cos(k*t+1);
fint = 2*sin(1)*sin(k)/k;
froots = ((-k:k)*pi-1)/k;
froots = froots(froots >= -1 & froots <= 1)';

N = 35;
L = legendre.matrix(N);
[x, w] = legendre.gauss(N);
c = L*f(x);

% Approximation
xtest = linspace(-1, 1, 1000)';
[P, D] = legendre.deriv_vec(N-1, xtest);
ferr = norm(P*c - f(xtest), inf) / norm(f(xtest), inf);
fperr = norm(D*c - fp(xtest), inf) / norm(fp(xtest), inf);
assert(ferr < 5e-14)
assert(fperr < 1e-12)

% Quadrature
Q = sum(f(x).*w);
quaderr = abs(Q-fint) / abs(fint);
assert(quaderr < 5e-15)

% Roots
r = legendre.roots(c);
r = r(abs(imag(r))<eps & real(r) >= -1 & real(r) <= 1);
r = sort(real(r));
rooterr = norm(r - froots, inf);
assert(rooterr < 5e-15)

