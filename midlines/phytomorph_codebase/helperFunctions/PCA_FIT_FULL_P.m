function [S C U E L ERR LAM] = PCA_FIT_FULL_P(M,COM)
%%%%%%%%%%%%%%%%
% INPUTS:   M       : = data matrix
%           COM     : = number of vectors 
%%%%%%%%%%%%%%%%
% OUTPUTS:  S       : = simulated signal (reconstruction)
%           C       : = components ("unique" fingerprint)
%           U       : = mean of data
%           E       : = basis vectors 
%           L       : = eig values
%           ERR     : = error in reconstruction
%           LAM     : = percent explained
%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % performs PCA (principal component analysis) using eigenvector
    % decomposition, backprojects to simulate the data and calculates the
    % error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:   
    %           M       : = data matrix
    %           COM     : = number of vectors 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           S       : = simulated signal (reconstruction)
    %           C       : = components ("unique" fingerprint)
    %           U       : = mean of data
    %           E       : = basis vectors 
    %           L       : = eig values
    %           ERR     : = error in reconstruction
    %           LAM     : = percent explained
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % take the mean
    fprintf(['PCA:start:taking off mean \n']);
    toOpDim = 1;
    U = mean(M,toOpDim);
    M = bsxfun(@minus,M,U);
    fprintf(['PCA:end:taking off mean \n']);
    % look at covariance
    fprintf(['PCA:start:creating COV \n']);
    
    COV = zeros(size(M,2));
    for e1 = 1:size(M,2)
        tmp = M(:,e1);
        parfor e2 = 1:size(M,2)
             COV(e1,e2) = mtimesx(tmp,'T',M(:,e2),'speed');
        end
    end
    
    try
        COV = mtimesx(M,'T',M,'speed');
    catch
        COV = M'*M;
    end
    COV = COV / size(M,toOpDim);
    fprintf(['PCA:end:creating COV \n']);
    % eig vector decomp
    fprintf(['PCA:start:decomposing COV \n'])
    [E L] = eigs(COV,COM);
    fprintf(['PCA:end:decomposing COV \n'])
    % return eig values
    LAM = L;
    % calc percent explained
    L = cumsum(diag(L))*sum(diag(L))^-1;
    % get coeffs (fingerprints)
    try
        C = mtimesx(M,E);
    catch
        C = M*E;
    end
    % use the back-projection to create simulated signal
    S = PCA_BKPROJ(C,E,U);
    % calc the error
    ERR = sum((S - M).^2,1).^.5;


%{
% take the mean
fprintf(['PCA:start:taking off mean \n']);
U = mean(M,1);
for i = 1:size(M,1)
    M(i,:) = M(i,:) - U;
end
fprintf(['PCA:end:taking off mean \n']);
% look at covariance
COV = cov(M);
% eig vector decomp
[E L] = eigs(COV,COM);
% return eig values
LAM = L;
% calc percent explained
L = cumsum(diag(L))*sum(diag(L))^-1;
% get coeffs (fingerprints)
if isdeployed
    C = M*E;
else
    C = mtimesx(M,E);
end
% use the back-projection to create simulated signal
S = PCA_BKPROJ(C,E,U);
% calc the error
ERR = sum((S - M).^2,2).^.5;
    %}