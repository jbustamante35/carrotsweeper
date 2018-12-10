classdef curveSegment < handle
    
    properties
        orientation;
        segment;
        image;
        Osegment;
        imageName;
        centerPoint;
        labelVec;
    end
    
    methods
        
        function [obj] = curveSegment()
            obj = obj;
        end
        
        
        function [D] = distance(obj,target)
            rot = obj.orientation(:,1)'*target.orientation(:,1);
            bend = mean(sum((obj.Osegment - target.Osegment).^2,1)).^.5;
            imgD = mean(abs(obj.image(:) - target.image(:)));
            D(1) = rot;
            D(2) = bend;
            D(3) = imgD;
        end
    end
    
end