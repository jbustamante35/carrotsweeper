classdef phytoPatch  < phytoGeo
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
        function [obj] = phytoPatch(varargin)
            % super constructor
            obj = obj@phytoGeo();            
            % set default views            
            obj.view_props.type = 'phytoPatch';
            % init point(s)
            if nargin == 1
                % init vars
                for e = 1:nargin
                    if isa(varargin{e},'phytoDomain')
                        % assign domain
                        obj.setDomain(varargin{e});
                    else
                        % assign data
                        S(1).type = '()';
                        S(1).subs{1} = 1:1:size(varargin{1},1);
                        S(1).subs{2} = ':';
                        subsasgn(obj,S,varargin{1});                        
                    end
                end
            end
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