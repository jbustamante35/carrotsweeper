classdef myTransitionFunction
    properties
        
    end
    
    methods (Abstract)
        [prob] = computeTransitionProb(obj,stateObservation,augVars);
        
        [] = updateTransitionData(obj,updatedData);
    end
end