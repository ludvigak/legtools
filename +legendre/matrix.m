function L = matrix(order, args)
    arguments
        order
        args.vpa = false
    end
    % L = legendre_matrix(order)
    % 
    % L computes Legendre expansion coefficients for a function defined at the
    % Gauss-Legendre quadrature nodes on [-1, 1]
    %
    %  c = L*f
    %
    % c_l = (2l+1)/2 \sum_n P_l(x_n) f_n w_n
    
    [x, w] = legendre.gauss(order, vpa=args.vpa);
    P = legendre.vec(order-1, x); %P = [P_0(x) P_1(x) P_2(x) ...]
    l = 0:order-1;
    %L = bsxfun(@times, P, w);
    L = P .* w;
    %L = bsxfun(@times, (2*l+1)/2, L);
    L = (2*l+1)/2 .*  L;
    L = L.';
end
