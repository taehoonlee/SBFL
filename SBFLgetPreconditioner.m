function E = SBFLgetPreconditioner(solvertype, l1, l2, row, col)
    switch solvertype
        case 'DCT'
            E = repmat(2*(1-cos(pi*(1:col)/(col+1))),row,1) + repmat(2*(1-cos(pi*(1:row)/(row+1)))',1,col);
            %E = ( sqrt(l1) + sqrt(l2) * E ) / row; % E = l1 * row * speye(row, col) + l2 * E; % miss of old
            % E = l1 * row * speye(row, col) + l2 * E; % miss of old
            E = ( l1 + l2 * E ) / row;
            %E = sqrt( ( l1 + l2 * E ) / row );
        case 'DFT'
            E = repmat(2*(1-cos(2*pi*(0:col-1)/col)),row,1) + repmat(2*(1-cos(2*pi*(0:row-1)/row))',1,col);
            %E = repmat(2*(1-cos(2*pi*(1:col)/col)),row,1) + repmat(2*(1-cos(2*pi*(1:row)/row))',1,col);
            % E = l1 * row * speye(row, col) + l2 * E; % miss of old
            E = ( l1 + l2 * E ) / row;
            %E = sqrt( ( l1 + l2 * E ) / row );
        case 'normalDST'
            E = 2 * diag(ones(1,row)) - diag(ones(1,row-1),1) - diag(ones(1,row-1),-1);
            E = kron(eye(col),E) + kron(E,eye(col));
            E = l1 * speye(row*col) + l2 * E;
        case 'normalDCT'
            E = 2 * diag(ones(1,row)) - diag(ones(1,row-1),1) - diag(ones(1,row-1),-1); E(1) = 1; E(end) = 1;
            E = kron(eye(col),E) + kron(E,eye(col));
            E = l1 * speye(row*col) + l2 * E;
        case 'normalDFT'
            E = 2 * diag(ones(1,row)) - diag(ones(1,row-1),1) - diag(ones(1,row-1),-1); E(row) = -1; E(end-row+1) = -1;
            E = kron(eye(col),E) + kron(E,eye(col));
            E = l1 * speye(row*col) + l2 * E;
        case 'chol'
            E = 2 * diag(ones(1,row)) - diag(ones(1,row-1),1) - diag(ones(1,row-1),-1); E(1) = 1; E(end) = 1;
            E = kron(eye(col),E) + kron(E,eye(col));
            E = l1 * eye(row*col) + l2 * E;
            E = chol(E, 'lower');
        otherwise
            E = [];
    end
end