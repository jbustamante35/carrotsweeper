function [curve sample] = extractCurves(img)
    % handle flip
    img = handleFLIP(img,[]);
    % refine contour
    B = contourRefine(I,1000);
end