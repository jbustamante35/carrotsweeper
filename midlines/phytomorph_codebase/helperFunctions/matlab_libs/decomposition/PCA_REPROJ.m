function [C] = PCA_REPROJ(M,E,U)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reproject 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:  
    %           M       : = data matrix
    %           E       : = basis vectors 
    %           U       : = mean of data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           C       : = unique finger print
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M = bsxfun(@minus,M,U);
    
    % subtract the mean
    %for i = 1:size(M,1)
    %    M(i,:) = M(i,:) - U;    % subtract the mean
    %end
    
    % project to the "smaller" - (rotate) - subspace
    C = M*E;                    % project to get the coeffs

end

