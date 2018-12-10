function [U,E,L] = PCA_FIT_FULLws(M,COM,sigma)
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
    % sigma default
    if nargin == 2
        sigma = 'largestabs';
    end
    % take the mean
    fprintf(['PCA:start:taking off mean \n']);
    toOpDim = 1;
    U = mean(M,toOpDim);
    M = bsxfun(@minus,M,U);
    fprintf(['PCA:end:taking off mean \n']);
    % look at covariance
    fprintf(['PCA:start:creating COV \n']);
    try
        fprintf(['i am speed.l.mcqueen \n']);
        COV = mtimesx(M,'T',M,'speed');
        fprintf(['and you know that.l.mcqueen \n']);
    catch
        COV = M'*M;
    end
    COV = COV / size(M,toOpDim);
    fprintf(['PCA:end:creating COV \n']);
    % eig vector decomp
    fprintf(['PCA:start:decomposing COV \n']);
    if ~isa(COV,'double')
        COV= double(COV);
    end
    [E,L] = eigs(COV,COM,sigma);
    fprintf(['PCA:end:decomposing COV \n']);
end
