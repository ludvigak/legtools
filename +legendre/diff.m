function d = diff(c)
% d = diff(c)
%
% Differentiate Legendre expansion

assert(size(c,2) == 1);

N = size(c, 1);
d = zeros(N-1, 1, 'like', c);

for i=1:N-1
    for j=i+1:2:N
        d(i) = d(i) + (2*i-1)*c(j);
    end
end