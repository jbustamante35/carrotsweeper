function [BV] = igetFrame(curve,radius)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculates the tangent and normal bundles for the curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:   
    %           curve   : = curve for analysis
    %           radius  : = radius for pca
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           BV      : = a sequence of basis vectors for the tangent and
    %           normal space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BV = [];
    for i = 1:size(curve,1)
        tBV = igetFrame_atP(curve,i,radius);
        BV = cat(3,BV,tBV);
    end
    % make the first dim the trial number
    BV = permute(BV,[3 1 2]);
end