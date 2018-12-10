classdef functionPredicate < functionObjectProps
    
    properties
        uniqueKey = '';
        OBJECT = 'RDF.PROPERTY';
        PREDICATE = 'RDFS.SUBPROPERTYOF';
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