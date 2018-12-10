function [D] = structureTensor(D,F)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from gradient information - create structure tensor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           D = stack of data - image with gradients belore
    %           F = filter parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           D = orginal information with the structure tensor glued
    %           onto bottom
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D1 = imfilter(D(:,:,1).^2,F);
    D2 = imfilter(D(:,:,2).^2,F);
    D3 = imfilter(D(:,:,1).*D(:,:,2),F);
    D = cat(3,D,D1,D2,D3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


