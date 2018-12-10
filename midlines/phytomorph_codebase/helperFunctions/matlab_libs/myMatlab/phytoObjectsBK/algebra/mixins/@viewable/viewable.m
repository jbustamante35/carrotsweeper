classdef viewable < matlab.mixin.Copyable & representation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        view_props = [];        % default view props
        view_rep = [];          % view rep
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [] = viewFree(obj,h,varargin)
            view(obj,h,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set view parameters
        function [] = setView(obj,para)
            obj.view_props = para;            
        end
    end
    methods (Abstract)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view
        [] = view(obj,h,frame,varargin);      
    end

    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parse view properties and attach to props
        function [props] = setProps(default,new)
            props = default;
            if ~isempty(new)
                flds = fields(new);
                for e = 1:numel(flds)
                    props.(flds{e}) = new.(flds{e});
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parse view properties and attach to props
        function [props] = parseView(v,props)
            for e = 1:numel(v)/2
                strt = (e-1)*2+1;
                ed = strt + 1;
                props.(v{strt}) = v{ed};
            end
        end
    end
    
end