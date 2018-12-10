function [S C U E L ERR LAM] = PCA_FIT_FULL_T(M,COM)
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
    toOpDim = 2;
    U = mean(M,toOpDim);
    M = bsxfun(@minus,M,U);
    fprintf(['PCA:end:taking off mean \n']);
    % look at covariance
    fprintf(['PCA:start:creating COV \n']);
    try
        COV = mtimesx(M,M,'T','speed');
    catch
        COV = M*M';
    end
    COV = COV / size(M,2);
    fprintf(['PCA:end:creating COV \n']);
    % eig vector decomp
    fprintf(['PCA:start:decomposing COV \n'])
    [E L] = eigs(COV,COM);
    [J sidx] = sort(diag(L),'descend');
    E = E(:,sidx);
    L = diag(L);
    L = L(sidx);
    L = diag(L);
    fprintf(['PCA:end:decomposing COV \n'])
    % return eig values
    LAM = L;
    % calc percent explained
    L = cumsum(diag(L))*sum(diag(L))^-1;
    % get coeffs (fingerprints)
    try
        C = mtimesx(E,'T',M);
    catch
        C = E'*M;
    end
    % use the back-projection to create simulated signal
    S = PCA_BKPROJ_T(C,E,U);
    % calc the error
    ERR = sum((S - M).^2,1).^.5;
end
