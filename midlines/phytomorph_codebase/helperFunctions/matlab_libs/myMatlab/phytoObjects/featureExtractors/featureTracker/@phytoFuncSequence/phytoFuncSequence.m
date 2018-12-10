classdef phytoFuncSequence  < matlab.mixin.Copyable
    
    properties
        para;
        notes;
        func;
    end
    
    methods
        %%%%%%%%%%%%
        % constructor
        function [obj] = phytoFuncSequence()
            obj.para = [];
            obj.notes = [];
            obj.func = [];
        end
        
        %%%%%%%%%%%%
        % put function in sequence
        function [para] = putFunc(obj,f,n)
            if nargin == 2;
                n = numel(obj.func);
            end
            obj.func{n} = f;
        end
        
        %%%%%%%%%%%%
        % get function from sequence
        function [ret] = getFunc(obj,n)
            ret = obj.func(n);
        end
        
        %%%%%%%%%%%%
        % eval chain
        function [X] = evalSequence(X)
            for f = 1:numel(obj.func)
                X = obj.func{f}(X);
            end
        end
        
    end
    
    methods (Static)
       
    end
    
    methods (Abstract)
       
    end
end