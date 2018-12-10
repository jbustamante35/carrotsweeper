function [angle,W,H] = getInformationSheetDimensions(points)

    hor1 = points(2,:,e) - points(1,:,e);
    hor2 = points(5,:,e) - points(4,:,e);
    ver1 = points(4,:,e) - points(1,:,e);
    ver2 = points(5,:,e) - points(2,:,e);

    W = mean([norm(hor1) norm(hor2)];
    H = mean([norm(ver1) norm(ver2)];
    
    
    angle = 180/pi*[atan2(-hor1(2),hor1(1)) atan2(-hor2(2),hor2(1)) atan2(-ver1(2),ver1(1))+pi/2 atan2(-ver2(2),ver2(1))+pi/2]
    angle = mean(angle);

end