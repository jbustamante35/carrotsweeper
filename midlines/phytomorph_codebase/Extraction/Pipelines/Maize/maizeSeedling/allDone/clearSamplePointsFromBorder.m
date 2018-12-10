function [pointList] = clearSamplePointsFromBorder(pointList,imageSZ,clearSZ,PD)
    rm = pointList(:,2) >= (imageSZ(1) - clearSZ(2)/2-PD);
    pointList(rm,:) = [];
    
    rm = pointList(:,1) >= (imageSZ(2) - clearSZ(1)/2-PD);
    pointList(rm,:) = [];
    
    rm = pointList(:,2) <= clearSZ(2)/2+PD;
    pointList(rm,:) = [];
    
    rm = pointList(:,1) <= clearSZ(1)/2+PD;
    pointList(rm,:) = [];
end