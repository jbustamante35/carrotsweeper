function [E L] = getSpan(M,COM)

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
    % look at covariance
    COV = M*M';
    % eig vector decomp
    [E L] = eigs(COV,COM);
end