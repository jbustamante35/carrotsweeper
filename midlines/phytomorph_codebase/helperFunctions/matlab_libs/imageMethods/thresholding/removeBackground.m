function [I BK Ibk] = removeBackground(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get background
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I       := image
    %           para    := parameters for running the script         
    %                   := para.sig         -> sigma for gaussian filter
    %                   := para.close       -> close amount
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           BK      := background       -> corner strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % get background
    BK = getBackground(I,para.getBackground);
    % remove background
    Ibk = I - BK;
end