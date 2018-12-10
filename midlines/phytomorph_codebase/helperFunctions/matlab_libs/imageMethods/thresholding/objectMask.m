function [B,thresh] = objectMask(I,thresh,closeVALUE,areaTHRESH)    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % objectMask - make mask for objects
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           I       : = image for interpolation
    %           para    : = parameters for make making
    %               para.filter.compute 
    %               para.filter.para
    %               para.BKsub.compute
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:
    %           Tout    = reference frame made up by the eigenvectors
    %           Y       = sampling in the reference frame Tout
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% makemask
    % ABOUT         := method for making mask
    % I             := image data
    % closeVALUE    := value for closing the image - for obtaining background 
    % areaTHRESH    := threshold value for image - if non - auto threshold
    if (size(I,3) ==3)
        I = rgb2gray(I);
    end
    
    %%%%%%%%%%
    % filter
    h = fspecial('gaussian',[10 10],4);
    I = imfilter(I,h,'replicate');
    
    %%%%%%%%%%
    % subtract the backgroud
    if ~isempty(closeVALUE)
        para.getBackground.close = closeVALUE;
        para.getBackground.sig = 2;
        I = removeBackground(I,para); 
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