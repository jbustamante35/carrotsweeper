classdef phytoApatch < myTb %< phytoGeo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % my node class
    properties
        domain;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoApatch(varargin)
            % super constructor
            obj = obj@myTb();            
            % set default views
            obj.view_props.props = [];
            obj.view_props.type = 'phytoApatch';
            % init point(s)
            if nargin == 1
                % init vars
                for e = 1:nargin
                    if isa(varargin{e},'phytoAdomain')
                        % assign domain
                        obj.setDomain(varargin{e});
                    else
                        obj.setData(varargin{e});                     
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,frame,vProps)
            S.type = '()';
            S.subs = {};            
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);            
            tensorView(subsref(obj,S),uProps,h);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set domain
        function [] = setDomain(obj,domain)
            obj.domain = domain;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get domain
        function [ret] = getDomain(obj)
            ret = obj.domain;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normalize
        function [] = normalize(obj,degree)
            
        end
    end
end