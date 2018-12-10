function [e] = errorMeasure(d1,d2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error measurement for two sets of vectors d1,d2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           d1      := first set of vectors
    %           d2      := second set of vectors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           e       := error measurement 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    e = (d1 - d2);
    e = sum(e.*e,2).^.5;
end