classdef sampleable < matlab.mixin.Copyable

    properties
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = sampleable(varargin)
            
        end
    end
    
    methods (Abstract)        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample object at domain
        [r] = sample(obj,image,domain);
    end
    
end