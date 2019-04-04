function r = roots(coeff)
% r = legendre.roots(coeff)
%
% Legendre expansion roots, sorted by ascending Bernstein radius

% Get roots
r = eig(legendre.comrade_matrix(coeff));

% Transform Berstein ellipses to horizontal lines
rmap = acos(r);

% Sort by imaginary component
[~, idx] = sort(abs(imag(rmap)));
r = r(idx);

