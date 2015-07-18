clear all;
close all;

l1 = 1;
l2 = 6;
mu1 = 1;
mu2 = 1;

% Data initialization: Lena
row = 64; col = 64; TOL = 1e-6; MAXITER = 200;
truebeta_mat = double(imresize(imread('Lena512.png'),[row col]));
truebeta = truebeta_mat(:);
y = truebeta + 8*randn(row*col, 1);

% Neighboring matrix initialization
L1 = diag(ones(1,row*col-1),1) - diag(ones(1,row*col)); L1(row:row:end,:) = [];
L3 = diag(ones(1,row*col-row),row) - diag(ones(1,row*col)); L3(end-row+1:end,:) = [];
L = sparse([L1; L3]);

beta1 = SBFL( [], y, l1, l2, mu1, mu2, row, col, L, TOL, MAXITER );
beta2 = SBFL_PCG( [], y, l1, l2, mu1, mu2, row, col, L, TOL, MAXITER, 'DFT' );
beta3 = SBFL_CGLS( [], y, l1, l2, mu1, mu2, row, col, L, TOL, MAXITER );
beta4 = SBFL_PCGLS( [], y, l1, l2, mu1, mu2, row, col, L, TOL, MAXITER, 'DFT' );

figure;
subplot(2,4,2),subimage(reshape(uint8(truebeta),row,col)); title('true beta'); axis off;
subplot(2,4,3),subimage(reshape(uint8(y),row,col)); title('y'); axis off;
subplot(2,4,5),subimage(reshape(uint8(beta1),row,col)); title('SBFL'); axis off;
subplot(2,4,6),subimage(reshape(uint8(beta2),row,col)); title('SBFL + PCG'); axis off;
subplot(2,4,7),subimage(reshape(uint8(beta3),row,col)); title('SBFL + CGLS'); axis off;
subplot(2,4,8),subimage(reshape(uint8(beta4),row,col)); title('SBFL + PCGLS'); axis off;