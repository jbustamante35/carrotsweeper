classdef sigmoidTransitionFunction < handle %< myTransitionFunction
    properties
        para = [];
        inv = 0;
    end
    
    methods
        
        function [obj] = sigmoidTransitionFunction(para,inv)
            obj.para = para;
            obj.inv = inv;
        end
        
        function [prob] = computeTransitionProb(obj,stateObservation,augVars)
            [~,prob] = mySigmoid_ver0(augVars{2},obj.para);
            if obj.inv
                prob = 1 - prob;
            end
        end
        
        function [] = updateTransitionData(obj,updatedData)
           
        end
    end
end