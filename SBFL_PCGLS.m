function [beta, f, residual, iter, f_inner] = SBFL_PCGLS(X, y, mu1, mu2, l1, l2, row, col, L, TOL, MAXITER, solvertype)

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
        A = [ eye(length(beta)); sqrt(l1) * eye(length(beta)); sqrt(l2) * L ];
    else
        A = [                 X; sqrt(l1) * eye(length(beta)); sqrt(l2) * L ];
    end

    for iter = 1:MAXITER

        if isempty(X)
            f(iter) = 1/2*norm(beta-y).^2 + l1*norm(beta,1) + l2*norm(L*beta,1);
            residual(iter) = norm(beta-y).^2;
        else
            f(iter) = 1/2*norm(X*beta-y).^2 + l1*norm(beta,1) + l2*norm(L*beta,1);
            residual(iter) = norm(X*beta-y).^2;
        end

        B = [y; 1 / sqrt(l1) * (l1*a - u); 1 / sqrt(l2) * (l2*b - v)];
        %beta = pcgls(A, eye(numel(beta)), B, zeros(size(beta)));
        
        x = zeros(size(beta));
        q = A' * B;
        z = SBFLlinsolve(q, row, col, solvertype, E);
        p = z;
        gamma0 = z' * z;
        gamma = gamma0;
        
        for k = 1:MAXITER

            if iter <= size(f_inner,1)
                f_inner(iter,k) = norm(gamma,2);
            end
            t = SBFLlinsolve(p, row, col, solvertype, E);
            Ap = A * t;
            alpha = gamma / ( Ap' * Ap );
            x = x + alpha * t;
            B = B - alpha * Ap;
            q = A' * B;
            z = SBFLlinsolve(q, row, col, solvertype, E);
            p = z + z'*z / gamma * p;
            gamma = z'*z;
            if norm(gamma,2) < 1e-6
                break;
            end

        end
        
        beta = x;

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