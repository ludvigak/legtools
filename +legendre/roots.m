function r = roots(coeff)
% Legendre expansion roots, sorted by ascending Berstein radius

% Get roots
r = eig(legendre.comrade_matrix(coeff));

% Transform Berstein ellipses to horizontal lines
rmap = acos(r);

% Sort by imaginary component
[~, idx] = sort(abs(imag(rmap)));
r = r(idx);

