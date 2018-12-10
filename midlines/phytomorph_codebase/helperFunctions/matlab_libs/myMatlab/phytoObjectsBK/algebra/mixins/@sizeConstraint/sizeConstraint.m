classdef sizeConstraint < constraint

    properties
        allowSize;
    end
    
    methods
        
        function [obj] = sizeConstraint(varargin)
            obj.allowSize = inf;
            if nargin == 1
                obj.allowSize = varargin{1};
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add elements to allowable list
        function [] = setSize(obj,sz)
            obj.allowSize = sz;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add elements to allowable list     
        function [b] = isAllow(obj,S)
            b = S(1).subs{1} <= obj.allowSize;
        end
        
        
    end
end
