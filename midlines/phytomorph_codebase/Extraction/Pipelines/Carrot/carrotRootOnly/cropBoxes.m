function [croppedBoxes] = cropBoxes(wholeImage,oPath,disp)
    
    % add a 50 pixel white border outide I of left and right
    I = padarray(wholeImage,[0 50],255);
    
    % Get enough of the black border around the boxes
    blackBoxMask = I(:,:,1) > 5 & I(:,:,1) < 105 & I(:,:,2) > 5 & I(:,:,2) < 105 & I(:,:,3) > 5 & I(:,:,3) < 105; % 100 for normal images, 105 for 9_rep1_14.47.27
    
    % Erode a little bit to seperate them from themselves and the border
    isolateBoxes = imerode(blackBoxMask, strel('square', 5));
    
    % Get rid of the border
    clearBorder = imclearborder(isolateBoxes);
   
    % Start filling in the borders (dilates and erode) 
    closeBoxes = imclose(clearBorder, strel('disk', 10, 0));
    
    % Finish filling in the borders (20 for normal images, 25 for
    % problematic 2-carrot pics)
    bridgeGaps = imdilate(closeBoxes, strel('line',20,0)); 
    
    % Not needed except for 8_rep1_14.41.35
    % bridgeGaps2 = imdilate(bridgeGaps, strel('line',5,90)); 
    
    % Fill in the boxes
    fillBoxes = imfill(bridgeGaps, 'holes');
    
    % Get rid of the residuals pixels
    cleanImage = bwareaopen(fillBoxes, 2700000);
    % imshow(cleanImage)

    R = regionprops(cleanImage, 'boundingbox');
    for b = 1:numel(R)
        croppedBoxes{b} = imcrop(I, R(b).BoundingBox);
    end
end