% Find complex roots of a simple function

addpath ..
    
f = @(t) 1 + 1i*t.^2;

N = 3;
[t,w] = legendre.gauss(N);
c = legendre.matrix(N)*f(t);

r = legendre.roots(c);
assert(norm(f(r), inf) < 5e-15);

zr = 0.2+1i;
tinit = 0.1;
tol = eps;
maxiter = 20;
[t, relres, iter] = legendre.newton(c, zr, tinit, maxiter, tol);
assert(abs(f(t) - zr) < 1e-14)