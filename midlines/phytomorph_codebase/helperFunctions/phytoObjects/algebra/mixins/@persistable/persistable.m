classdef persistable < matlab.mixin.Copyable
    
    properties 
        oStore;
    end
    
    methods (Abstract)        
        persist(obj);
        toBson(obj,xferO);        
    end
    
    methods (Static)        
        fromBson(varargin);        
    end
end