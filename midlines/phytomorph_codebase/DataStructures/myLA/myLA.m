classdef myLA < matlab.mixin.Heterogeneous & handle
   properties
       byte;
   end
   
   methods
        function [obj] = myLA(ibyte)
            if nargin == 0
                ibyte = 0;
            end
            obj.byte = logical(ibyte);
        end
        
        function [r] = eq(objA,objB)
            r = all(objA.byte == ~objB.byte);
        end
        
        function [] = glue(objA,objB)
            objA.byte = [objA.byte objB.byte];
        end
        
        function [] = swap(obj,p)
            for e = 1:numel(obj)
                obj(e).byte = obj(e).byte(p);
            end
        end
        
        
        function [] = mtimes(objA,objB)
            here =1;
        end
   end
end

%{
   
    q1 = myLA([1]);
    q2 = myLA([1]);
    qT = [q1 q2];


%}