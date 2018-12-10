classdef heavisideTransitionFunction < myTransitionFunction
    
    properties
        stateCount = [];
        switchFunction;
    end
    
    methods
        function [obj] = heavisideTransitionFunction(stateCount,switchFunction)
            obj.stateCount = stateCount;
            obj.switchFunction = switchFunction;
        end
        
        function [prob] = computeTransitionProb(obj,stateObservation,augVars)
            if obj.switchFunction(sum(stateObservation == augVars{1}),obj.stateCount)
                prob = 1;
            else
                prob = 0;
            end
        end
        
        function [] = updateTransitionData(obj,updatedData)
            
        end
    end
end