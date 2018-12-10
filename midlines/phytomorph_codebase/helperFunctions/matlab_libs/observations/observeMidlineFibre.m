function [L T] = observeMidlineFibre(M)
    %%%%%
    % derivative of midline
    dL = diff(M,1,1);
    dL = sum(dL.*dL,2).^.5;
    L = sum(dL);
    
    %%%%%
    % construct frame
    T = igetFrame_atP(M,1,30);
end