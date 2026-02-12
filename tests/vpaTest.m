% Repetition of some things from approxTest
addpath ..
clear

digits(32)

N = 40;
L = legendre.matrix(N,     vpa=true);
[x, w] = legendre.gauss(N, vpa=true);

xs = (x+1)/2;
ws = w/2;

gaussian = @(x) 2/sqrt(pi)*exp(-x.^2);
a = L*gaussian(xs);

resolution = double(norm(a(end-1:end), inf))
assert(resolution < 1e-25)

b = legendre.cumsum(a);
cumsum_diff_err = norm(legendre.diff(b)-a, inf) / norm(a, inf)
assert(cumsum_diff_err < 1e-32)


% Compute erf(1) - integral from 0 to 1
int_err = legendre.vec(N, 1)*b - ws' * legendre.vec(N-1, x)*a
assert(int_err < 1e-32)

