classdef closedCurveSegment < curveSegment
    properties
        centerVec;
        centerDistance;
        
        skeletonVec;
        skeletonDistance;
    end
    
    methods
        function [obj] = closedCurveSegment()
            obj = obj@curveSegment();
        end
    end
end