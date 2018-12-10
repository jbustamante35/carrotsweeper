function [qrA,qrW,qrH,whA,whW,whH,huA,huW,huH] = getInformationSheetDimensions(points)
    for e = 1:size(points,3)
        hor1 = points(2,:,e) - points(1,:,e);
        hor2 = points(5,:,e) - points(4,:,e);
        ver1 = points(4,:,e) - points(1,:,e);
        ver2 = points(5,:,e) - points(2,:,e);

        qrW(e) = mean([norm(hor1) norm(hor2)]);
        qrH(e) = mean([norm(ver1) norm(ver2)]);

        angle = 180/pi*[atan2(-hor1(2),hor1(1)) atan2(-hor2(2),hor2(1)) atan2(-ver1(2),ver1(1))+pi/2 atan2(-ver2(2),ver2(1))+pi/2];
        qrA(e) = mean(angle);
        
        hor1 = points(3,:,e) - points(1,:,e);
        hor2 = points(8,:,e) - points(7,:,e);
        ver1 = points(7,:,e) - points(1,:,e);
        ver2 = points(8,:,e) - points(3,:,e);

        whW(e) = mean([norm(hor1) norm(hor2)]);
        whH(e) = mean([norm(ver1) norm(ver2)]);
        
        angle = 180/pi*[atan2(-hor1(2),hor1(1)) atan2(-hor2(2),hor2(1)) atan2(-ver1(2),ver1(1))+pi/2 atan2(-ver2(2),ver2(1))+pi/2];
        whA(e) = mean(angle);
        
        
        hor1 = points(3,:,e) - points(2,:,e);
        hor2 = points(6,:,e) - points(5,:,e);
        ver1 = points(5,:,e) - points(2,:,e);
        ver2 = points(6,:,e) - points(3,:,e);

        huW(e) = mean([norm(hor1) norm(hor2)]);
        huH(e) = mean([norm(ver1) norm(ver2)]);
        
        angle = 180/pi*[atan2(-hor1(2),hor1(1)) atan2(-hor2(2),hor2(1)) atan2(-ver1(2),ver1(1))+pi/2 atan2(-ver2(2),ver2(1))+pi/2];
        huA(e) = mean(angle);

       
    end

end