function b = cumsum(a)
% b = coefficients of F(x) = \int_0^x f(t) dt
% f(x) = sum a_n P_n(x)
%
% a, b are (N+1)xK sym matrices allowed

    % C = integration constant
    C = 0;

    N = size(a,1) - 1;
    K = size(a,2);

    b = zeros(N+2, K, class(a));   % degree increases by 1

    % Constant of integration
    b(1,:) = C;

    % Contribution from a_0: ∫P_0 = x = P_1
    if N >= 0
        b(2,:) = b(2,:) + a(1,:);
    end

    % n >= 1 terms
    for n = 1:N
        coeff = a(n+1,:) / (2*n+1);
        b(n,:)   = b(n,:)   - coeff;   % P_{n-1}
        b(n+2,:) = b(n+2,:) + coeff;   % P_{n+1}
    end
end
