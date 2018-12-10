classdef variable < functionObjects
    
    properties
        
    end
    
    methods
        function [obj] = variable(value,valueDatabase,valueStream,toPersist)
            obj = obj@functionObjects();
        end
    end
    
    methods (Access = private)
    end
end