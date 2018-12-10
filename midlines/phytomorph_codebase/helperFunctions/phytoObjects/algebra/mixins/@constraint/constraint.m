classdef constraint < matlab.mixin.Copyable

    % my tensor class
    properties
       
    end
    
    methods (Abstract)
        [b] = isAllow(obj,ele);
    end

end