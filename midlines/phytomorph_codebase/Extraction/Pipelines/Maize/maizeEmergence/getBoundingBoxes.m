function [R] = getBoundingBoxes(MASK)
    R = regionprops(logical(MASK),'BoundingBox','Centroid');
end