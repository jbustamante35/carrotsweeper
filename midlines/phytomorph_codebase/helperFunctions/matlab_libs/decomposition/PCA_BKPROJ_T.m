function M = PCA_BKPROJ_T(C,E,U)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % back project
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:   
    %           C       : = components ("unique" fingerprint)
    %           E       : = basis vectors 
    %           U       : = mean of data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           M       : = create simulated data - from model -> create "data"
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % coeffs
    %M = (E*C);                            % project to get the coeffs
    try
        M = mtimesx(E,C);
    catch
        M = E*C;
    end
    M = bsxfun(@plus,M,U);
end
