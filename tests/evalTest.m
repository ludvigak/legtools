% Test Legendre poly evaluation across implementations
addpath ..
    
rng(1);
n = 10;
N = 20;
x = (1-2*rand(n,1)) + 1i*(1-2*rand(n,1));

% Vectorized in argument
P1 = legendre.vec(N, x);
[P2, D2] = legendre.deriv_vec(N, x);
% Scalar argument
P3 = zeros(size(P1));
P4 = zeros(size(P1));
D4 = zeros(size(P1));
for i=1:n
    P3(i,:) = legendre.scalar(N, x(i));
    [P4(i,:), D4(i,:)] = legendre.deriv_scalar(N, x(i));
end
% Compare
assert(norm(P1-P2,inf) < eps)
assert(norm(P1-P3,inf) < eps)
assert(norm(P1-P4,inf) < eps)
assert(norm(D2-D4,inf) < eps)

% Compare to native Matlab 
x = linspace(-1,1, n);
P = legendre.vec(N, x);
for l=0:N
    p = legendre(l, x);
    assert(norm(p(1,:)' - P(:,l+1), inf) < 2*eps)
end
