function [C ERR mu sd] = PCA_REPROJ_T(M,E,U)
    ERR = [];
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
    fprintf(['PCA_RE:start:taking off mean \n']);
    M = bsxfun(@minus,M,U);
    fprintf(['PCA_RE:end:taking off mean \n']);

    % subtract the mean
    %for i = 1:size(M,1)
    %    M(i,:) = M(i,:) - U;    % subtract the mean
    %end
    fprintf(['PCA_RE:start:proj into \n']);
    % project to the "smaller" - (rotate) - subspace
    try
        C = mtimesx(E,'T',M);                    % project to get the coeffs
    catch
        C = E'*M;
    end
    fprintf(['PCA_RE:end:proj into \n']);
    if nargout >= 2
        Mi = PCA_BKPROJ_T(C,E,U);
        ERR = sum((bsxfun(@plus,M,U) - Mi).^2,1).^.5;
    end
    if nargout > 2
        delta = (M - Mi);
        mu = mean(delta,1);
        sd = std(delta,1,1);
    end
end

