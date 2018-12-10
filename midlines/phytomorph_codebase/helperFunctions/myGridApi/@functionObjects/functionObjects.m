classdef functionObjects < functionObjectProps
    
    
    properties
        OBJECT = 'RDFS.CLASS';
        PREDICATE = 'RDFS.SUBCLASSOF';
    end
    
    methods
        function [obj]= functionObjects()
           obj = obj@functionObjectProps();
        end
        
        
    end
    
    methods 
        
    end
    
    methods (Static)
        function [classDefURL] = renderClassDef(obj,valueDatabase)
            classDefURL = functionObjectProps.renderClassDef(obj,valueDatabase,functionObjects.OBJECT);
        end
        
        function [] = renderClassDefLevels(obj,valueDatabase)
            functionObjectProps.renderClassDefLevels(obj,valueDatabase,functionObjects.OBJECT,functionObjects.PREDICATE);
        end
    end
end