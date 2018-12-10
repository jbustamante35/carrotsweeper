function [MSK] = makeBoundaryMask(sz,mskValue)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make Boundary Mask Image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           sz          := size of mask
    %           mskValue    := border values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           MSK = mask, to exclude the rim of the images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make mask to exclude points within radius of the image boundary. 
    MSK = ones(sz);
    % left
    if mskValue(1) > 1 ; MSK(1:mskValue(1),:)       = 0; end    
    % right
    if mskValue(2) > 1 ; MSK(end-mskValue(2):end,:) = 0; end
    % top
    if mskValue(3) > 1 ; MSK(:,1:mskValue(3))       = 0; end
    % bottom
    if mskValue(4) > 1 ; MSK(:,end-mskValue(4):end) = 0; end    
end
    
