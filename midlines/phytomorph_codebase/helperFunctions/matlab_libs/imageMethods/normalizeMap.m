function [nI para] = normalizeMap(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalize the image - if para is given then use it to normalize
    % if para is not given then min,max to return para for next input.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I       := vector to normalize
    %           para    := parameters for normalizing 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           nI      := normalized image
    %           para    := parameters for normalizing 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (nargin == 1) || isempty(para)
        para(1) = min(I(:));
        nI = I - para(1);
        para(2) = max(nI(:))^-1;
        nI = nI * para(2);
    else        
        nI = I - para(1);
        nI = nI * para(2);        
    end
end