classdef geometrical < matlab.mixin.Copyable
    % each geometry object hasa tensor which
    % is second rank and maps between the natural frame and
    % another frame.  The natural frame or object frame is
    % a property of the object.
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%
        bf;             % basis frame
        %%%%%%%%%%%%%%%%%%%%%%%%
        s;              % state config for default decompose
                        % if defined then prefer else not yes
        %%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = geometrical(varargin)
            obj.s = 1;
            obj.bf = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set normalization state
        function [] = setState(obj,s)
            obj.s = s;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set basis functions
        function [] = setBasis(obj,b)
            obj.bf = b;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get basis functions
        function [b] = getBasis(obj)
            b = obj.bf;
        end
    end
end