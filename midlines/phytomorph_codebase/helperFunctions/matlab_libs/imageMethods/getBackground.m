function [BK] = getBackground(I,para)
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
    %%% filter
    h = fspecial('gaussian',[3*para.sig 3*para.sig],para.sig);    
    I = imfilter(I,h,'replicate');
    %%% close for background
    BK = imclose(I,strel('disk',para.close));
end