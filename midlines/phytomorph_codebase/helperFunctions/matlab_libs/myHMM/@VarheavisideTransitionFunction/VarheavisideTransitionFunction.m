classdef VarheavisideTransitionFunction < myTransitionFunction
    
    properties
        stateCount = [];
        switchFunction;
        outProb = 1;
        nodeIDX = 1;
    end
    
    methods
        function [obj] = VarheavisideTransitionFunction(stateCount,switchFunction,outProb,nodeIDX)
            obj.stateCount = stateCount;
            obj.switchFunction = switchFunction;
            obj.outProb = outProb;
            obj.nodeIDX = nodeIDX;
        end
        
        function [prob] = computeTransitionProb(obj,stateObservation,augVars)
            if numel(stateObservation) > 0
                st = stateObservation(end);
                R = regionprops(stateObservation(end) == stateObservation,'PixelIdxList');
                stateObservation = zeros(size(stateObservation));
                stateObservation(R(end).PixelIdxList) = st;
            end
            if obj.switchFunction(sum(stateObservation == augVars{1}),obj.stateCount)
                prob = obj.outProb;
            else
                prob = 1-obj.outProb;
            end
        end
        
        function [] = updateTransitionData(obj,updatedData)
            
        end
    end
end