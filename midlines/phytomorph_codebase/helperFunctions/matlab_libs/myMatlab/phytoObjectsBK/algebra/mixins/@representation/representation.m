classdef representation < matlab.mixin.Copyable

    properties
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = representation(varargin)
            
        end
    end
    
    methods (Abstract)        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation
        [r] = rep(obj,frame,type);
    end
    
end