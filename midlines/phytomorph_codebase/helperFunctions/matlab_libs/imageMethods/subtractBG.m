function [I] = subtractBG(I)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subtract background of the respective image - 
    % background determined by a gaussian filtering followed by dilating
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT:    I  = image (unchanged)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           I  = image (background subtracted)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    h = fspecial('gaussian',7,3);
    I = imfilter(I,h);
    se = strel('disk',31);
    BK = imdilate(I,se);
    I = I - BK;
    I = normalizeImage(I);
end

