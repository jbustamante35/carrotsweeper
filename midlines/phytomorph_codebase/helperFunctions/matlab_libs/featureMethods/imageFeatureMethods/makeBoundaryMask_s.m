function [MSK] = makeBoundaryMask_s(sz,mskValue)    
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
    % top
    if mskValue(1) > 1 ; MSK(1:mskValue(1),:)       = 0; end    
    % middle
    if mskValue(2) > 1 ; MSK(:,((round(end/2))- mskValue(2)):((round(end/2))+ mskValue(2))) = 0; end
    % bottom
    if mskValue(3) > 1 ; MSK(end-mskValue(3):end,:) = 0; end
    % left
    if mskValue(4) > 1 ; MSK(:,1:mskValue(4))       = 0; end
    %right
    if mskValue(5) > 1 ; MSK(:,end-mskValue(5):end) = 0; end    
end
    
