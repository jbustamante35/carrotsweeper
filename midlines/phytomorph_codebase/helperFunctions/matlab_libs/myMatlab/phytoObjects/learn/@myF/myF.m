classdef myF < matlab.mixin.Copyable
    % my functional space
    properties    
        f;      % function list of type myiF
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myF(varargin)
            f = {};
        end        
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add functional to set
        function [obj] = addFunctional(obj,f,idx)
           if nargin <= 2
                idx = numel(obj.f) + 1;
            end
            obj.f{idx} = f;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % eval functional(s)
        function [ret] = e(obj,d)
            for fi = 1:numel(obj.f)
                ret{fi} = obj.f{fi}.e(d);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % invert functional(s)
        function [] = i(obj,d,c)
            for fi = 1:numel(obj.f)
                obj.f{fi}.i(d,c);
            end
        end
        
    end
end