classdef germCell
    
    properties
        chr = [];
        xm = [];
    end
    
    methods
        function [obj] = germCell()
            obj.chr = CH();
            obj.xm = xM(obj.chr);
        end
    end
end