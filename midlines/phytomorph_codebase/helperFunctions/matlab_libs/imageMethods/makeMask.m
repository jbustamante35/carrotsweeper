function [B,thresh] = makeMask(I,thresh,closeVALUE,areaTHRESH)
    % ABOUT         := method for making mask
    % I             := image data
    % closeVALUE    := value for closing the image - for obtaining background 
    % areaTHRESH    := threshold value for image - if non - auto threshold
    if (size(I,3) ==3)
        I = rgb2gray(I);
    end
    %%%%%%%%%%
    % subtract the backgroud
    if ~isempty(closeVALUE)
        %%% obtain the background
        BK = getBackground(I,closeVALUE);
        %%% subtract the background
        I = I - BK;
    end
    %%%%%%%%%%
    % normalize the background
    %%% normalize the image
    I = normalize(I);
    %%%%%%%%%%
    % obtain the threshold value    
    if isempty(thresh)
        thresh = graythresh(255*I);
    end
    %%%%%%%%%%
    % threshold    
    B = I < thresh;
    %%%%%%%%%%
    % auto fill holes    
    B = imfill(B,'holes');
    %%%%%%%%%%
    %%% remove small objects
    B = bwareaopen(B,areaTHRESH);
    %%%%%%%%%%
    %%% errode
    %B = imerode(B,strel('disk',5));
end