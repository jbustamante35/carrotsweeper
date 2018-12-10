function [coneTainer,coneTainerBB] = processMaskSet(I,M,bottomBuffer)
    coneTainer = M(:,:,2);
    coneTainer = bwareaopen(coneTainer,1000);
    coneTainer = imclose(coneTainer,strel('square',21));
    coneTainer = imfill(coneTainer,'holes');
    coneTainer(end-bottomBuffer:end,:) = 0;
    
    
    
    coneTainer = imopen(coneTainer,strel('square',31));
    
    
    coneTainer = bwlarge(coneTainer,3);
    
    coneTainerBB = regionprops(coneTainer,'BoundingBox','Centroid');
    
    
    
    
end