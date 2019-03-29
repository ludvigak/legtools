function A = comrade_matrix(coeff)
% Compute comrade matrix for polynomial given by Legendre
% expansion p(x) = \sum_k coeff_k*P_k(x)
% The roots of p(x) are the eigenvalues of A.

    N = numel(coeff)-1;
    [alpha, gamma] = recurrence_coeffs(N);
    A = zeros(N,N);
    %A(1,1) = beta(1);
    A(2,1) = alpha(1);
    for i=2:N-1
        A(i-1,i) = gamma(i);
        %A(i  ,i) = beta(i);
        A(i+1,i) = alpha(i);     
    end
    A(N-1,N) = gamma(N);
    %A(N  ,N) = beta(N);
    A(:,N) = A(:,N) - coeff(1:N) .* alpha(N)/coeff(N+1);
end
    
function [alpha, gamma] = recurrence_coeffs(N)
% Coefficients of three-term recurrence 
% alpha_k*P_{k+1} = (x-beta_k)*P_{k} - gamma_k*P_{k-2}
    k = 0:N-1;
    alpha = (k+1)./(2*k+1);
    % (beta=0)
    gamma = k./(2*k+1);
end