function [beta, f, residual, iter, f_inner] = SBFL(X, y, mu1, mu2, l1, l2, row, col, L, TOL, MAXITER)

    beta = zeros(size(L,2), 1);
    a = zeros(size(L,2), 1);
    b = zeros(size(L,1), 1);
    u = zeros(size(L,2), 1);
    v = zeros(size(L,1), 1);
    
    f = zeros(1, MAXITER);
    residual = zeros(1, MAXITER);
    f_inner = zeros(MAXITER, MAXITER);
    
    if isempty(X)
        A = speye(length(beta)) + l1 * speye(length(beta)) + l2 * (L' * L);
        XTy = y;
    else
        A = X' * X + l1 * speye(length(beta)) + l2 * (L' * L);
        XTy = X' * y;
    end

    for iter = 1:MAXITER

        if isempty(X)
            f(iter) = 1/2*norm(beta-y).^2 + l1*norm(beta,1) + l2*norm(L*beta,1);
            residual(iter) = norm(beta-y).^2;
        else
            f(iter) = 1/2*norm(X*beta-y).^2 + l1*norm(beta,1) + l2*norm(L*beta,1);
            residual(iter) = norm(X*beta-y).^2;
        end

        B = XTy + l1 * (a - u) + l2 * L' * (b - v);
        
        beta = A \ B;

        a = SBFLsoftThresholding(beta + 1 / l1 * u, mu1/l1);
        u = u + l1 * (beta - a);

        b = SBFLsoftThresholding(L * beta + 1 / l2 * v, mu2/l2);
        v = v + l2 * (L * beta - b);
        
        if iter > 1
            if abs(f(iter) - f(iter-1)) / f(iter-1) < TOL
                break;
            end
            %fprintf('%d %g\n',iter, abs(f(iter) - f(iter-1)) / f(iter-1));
        end

    end

end