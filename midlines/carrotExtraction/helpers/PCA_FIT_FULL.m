function [S, C, U, E, L, ERR, LAM] = PCA_FIT_FULL(M, COM)
%% PCA_FIT_FULL: Principal Components Analysis using eigenvector decomposition
% Performs PCA (principal component analysis) using eigenvector decomposition,
% then backprojects to simulate the data and calculates error.
%
% Usage: 
%   [S, C, U, E, L, ERR, LAM] = PCA_FIT_FULL(M, COM)
%
% Input:
%   M: data matrix
%   COM: number of vectors
%
% Output:
%   S: simulated signal (reconstruction)
%   C: components ("unique" fingerprint)
%   U: mean of data
%   E: basis vectors
%   L: eig values
%   ERR: error in reconstruction
%   LAM: percent explained
%

%% Get the mean
toOpDim = 1;
U       = mean(M, toOpDim);
M       = bsxfun(@minus,M,U);

%% Get covariance
try
    COV = mtimesx(M,'T',M,'speed');
catch
    COV = M' * M;
end

COV = COV / size(M, toOpDim);

%% Eigenvector decomposition
[E, L]    = eigs(COV, COM);
[~, sIdx] = sort(diag(L), 'descend');

%% Return eig values
E   = E(:, sIdx);
L   = diag(L);
L   = L(sIdx);
L   = diag(L);
LAM = L;

% Calculate percent explained
L = cumsum(diag(L)) * sum(diag(L))^-1;

%% Coeffs (fingerprints)
try
    C = mtimesx(M, E);
catch
    C = M * E;
end

%% Back-projection to create simulated signal and compute error
S   = PCA_BKPROJ(C, E, U);
ERR = sum((S - M).^2, 1).^0.5;

end