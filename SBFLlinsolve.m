function z = SBFLlinsolve(r, row, col, solvertype, L)
    switch solvertype
        case 'DCT'
            zz = fft2([ zeros(1,2*col+2); ...
                zeros(row,1), reshape(r,row,col), zeros(row,1), -fliplr(reshape(r,row,col)); ...
                zeros(1,2*col+2);
                zeros(row,1), -flipud(reshape(r,row,col)), zeros(row,1), rot90(reshape(r,row,col),2) ]);
            U = -real(zz(2:(row+1),2:(col+1))) / 4;
            Z = U./L';
            vv = fft2([ zeros(1,2*col+2); ...
                zeros(row,1), Z, zeros(row,1), -fliplr(Z); ...
                zeros(1,2*col+2); ...
                zeros(row,1), -flipud(Z), zeros(row,1), rot90(Z,2) ]);
            z = -real(vv(2:(row+1),2:(col+1))) / (row+1) / (col+1);
            z = reshape(z, [], 1);
        case 'DFT'
            U = fft2(reshape(r,row,col)) / sqrt(row*col); U(2:end,:) = U(end:-1:2,:); U(:,2:end) = U(:,end:-1:2);
            Z = U./L';
            z = real(fft2(Z)) / (row*col); % sqrt(row*col); % miss of old
            z = reshape(z, [], 1);
        case 'id'
            z = r / sqrt(row);
        otherwise
            z = L \ r;
    end
end