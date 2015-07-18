function [beta, f, residual, iter, f_inner] = SBFL_PCG(X, y, mu1, mu2, l1, l2, row, col, L, TOL, MAXITER, solvertype)

    beta = zeros(size(L,2), 1);
    a = zeros(size(L,2), 1);
    b = zeros(size(L,1), 1);
    u = zeros(size(L,2), 1);
    v = zeros(size(L,1), 1);
    
    f = zeros(1, MAXITER);
    residual = zeros(1, MAXITER);
    f_inner = zeros(MAXITER, MAXITER);
    E = SBFLgetPreconditioner(solvertype, l1, l2, row, col);
    
    if isempty(X)
        A = speye(length(beta)) + l1 * speye(length(beta)) + l2 * (L' * L);
        XTy = y;
    else
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

        if isempty(X)
            r = B - A * beta;
        else
            r = B - X' * (X * beta) - l1 * beta - l2 * (L' * L) * beta;
        end
        z = SBFLlinsolve(r, row, col, solvertype, E);
        p = z;
        rzold = r' * z;
        
        for k = 1:MAXITER

            if iter <= size(f_inner,1)
                f_inner(iter,k) = norm(r,2);
            end
            
            if isempty(X)
                Ap = A * p;
            else
                Ap = X' * (X * p) + l1 * p + l2 * (L' * L) * p;
            end
            alpha = rzold / (p'*Ap);
            beta = beta + alpha*p;
            r = r - alpha * Ap;
            z = SBFLlinsolve(r, row, col, solvertype, E);
            rznew = r'*z;
            p = z + rznew / rzold * p;
            rzold = rznew;
            if norm(r,2) < 1e-8
                break;
            end

        end

        a = SBFLsoftThresholding(beta + 1 / l1 * u, mu1/l1);%a = u + beta;
        u = u + l1 * (beta - a);

        b = SBFLsoftThresholding(L * beta + 1 / l2 * v, mu2/l2);%b = v + L * beta;
        v = v + l2 * (L * beta - b);
        
        if iter > 1
            if abs(f(iter) - f(iter-1)) / f(iter-1) < TOL
                break;
            end
            %fprintf('%d %g\n',iter, abs(f(iter) - f(iter-1)) / f(iter-1));
        end

    end

end