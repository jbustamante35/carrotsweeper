function [a] = measurePointBox(box)

    topHLine = box(2,:) - box(1,:);
    bottomHLine = box(4,:) - box(3,:);
    leftVLine = box(1,:) - box(3,:);
    rightVLine = box(2,:) - box(4,:);
    a(1) = atan2(topHLine(1),topHLine(2));
    a(2) = atan2(bottomHLine(1),bottomHLine(2));
    a(3) = atan2(leftVLine(1),leftVLine(2))+pi/2;
    a(4) = atan2(rightVLine(1),rightVLine(2))+pi/2;
    
    a = mean(a);
    
    
end