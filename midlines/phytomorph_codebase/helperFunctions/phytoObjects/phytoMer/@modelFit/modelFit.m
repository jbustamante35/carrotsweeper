classdef modelFit < matlab.mixin.Copyable
    properties
        parameters;
        Xdata;
        Ydata;
    end
    
    methods (Abstract)
        fit;
        
        fnval(obj,X);
        
    end
end