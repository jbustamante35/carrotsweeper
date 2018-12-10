function [baseSegment contourBody] = separateContour(contour)
    fidx = find(contour(:,2) == 1);
    
    baseSegment = contour(fidx,:);
    contourBody = contour(~fidx,:);
    
end