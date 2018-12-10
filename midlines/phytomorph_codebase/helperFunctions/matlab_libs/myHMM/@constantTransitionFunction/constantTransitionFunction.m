classdef constantTransitionFunction < handle %< myTransitionFunction
    properties
        transitionProb = [];
    end
    
    methods
        
        function [obj] = constantTransitionFunction(transitionProb)
            obj.transitionProb = transitionProb;
        end
        
        function [prob] = computeTransitionProb(obj,stateObservation,augVars)
            prob = obj.transitionProb;
        end
        
        function [] = updateTransitionData(obj,updatedData)
            obj.transitionProb = updatedData;
        end
    end
end