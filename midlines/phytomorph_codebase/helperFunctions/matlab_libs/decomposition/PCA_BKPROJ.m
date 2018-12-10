function M = PCA_BKPROJ(C,E,U)
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
    M = (E*C')';                            % project to get the coeffs
    M = bsxfun(@plus,M,U);
    %{
    for i = 1:size(M,1)
        M(i,:) = M(i,:) + U;                % subtract the mean
    end
    %}
end
