function [midline,contour] = getMidlineAndContour(carrotMask,ORIN)
    midline = [];
    contour = [];
    [out] = isolate_carrot_Roots(double(~carrotMask),0,[],[]);
    midline = out(1).midlines.data';
    sz = size(carrotMask);


    
    
    rm = midline(:,1) < 300;
    midline(rm,:) = [];
    midline(:,1) = midline(:,1) - 300;
    if ~ORIN
        midline(:,1) = midline(:,1) - sz(2)/2;
        midline(:,1) = -midline(:,1);
        midline(:,1) = midline(:,1) + sz(2)/2;
    end
    
    contour = out(1).contours.data';
    rm = contour(:,1) < 300;
    contour(rm,:) = [];
    contour(:,1) = contour(:,1) - 300;
    
    if ~ORIN
        contour(:,1) = contour(:,1) - sz(2)/2;
        contour(:,1) = -contour(:,1);
        contour(:,1) = contour(:,1) + sz(2)/2;
    end
    
end
