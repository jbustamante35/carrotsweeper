classdef myHS < matlab.mixin.Copyable & viewable
    % abstract def of hetero-geneous container
    % the subclasses have added constraints: types are allowed or not,
    % allowed types can include the type of the container. length or size
    % constraints can be added. emission constraints added.
    % of types: if there are two types, and one allowed type is the
    % container. then the other can be the element type. when
    % this happens, we can view there to be only one allowed ype, the
    % container. then the existance of a single element next to a set ina
    % set turns into two sets ina set with one containing one and one many
    % {e,{e,e,e}} --> {{e},{e,e,e}} the difference might be in the return
    % type s{1} is e vs {e}. this changes extraction to be size based and
    % typed based to ensure rather then type based only, 
    properties
        S = {};                 % set                     
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsref
        function [r] = subsref(obj,S)
            if strcmp(S(1).type,'.')
                % put note when this happens: when nargout == 0 and error
                % re return value
                if nargout == 0
                    % ans: when there is a function call and S is numel==2
                    if numel(S) == 1
                        ans =  builtin('subsref',obj,S);
                        ans
                    else
                        builtin('subsref',obj,S);
                    end
                else
                    r = builtin('subsref',obj,S);
                end
            elseif strcmp(S(1).type,'()')
                r = builtin('subsref',obj.S,S); 
            elseif strcmp(S(1).type,'{}')
                e = builtin('subsref',obj.S,S);
                r = e;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [r] = subsasgn(obj,S,B)
            if strcmp(S(1).type,'.')
                r = builtin('subsasgn',obj,S,B);
            elseif strcmp(S(1).type,'()')
                obj.S = builtin('subsasgn',obj.S,S,B);
                r = obj;
            elseif strcmp(S(1).type,'{}')
                obj.S = builtin('subsasgn',obj.S,S,B);                
                r = obj;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overloaded numel
        function [ret] = numel(obj,varargin)
            if nargin == 1
                ret = numel(obj.S);
            else
                ret = builtin('numel',obj,varargin{1});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % type of elements in set - first level deep
        function [ret] = type(obj)
            ret = {};
            S(1).type = '{}';
            S(1).subs{1} = 1;
            for e = 1:numel(obj)
                ele = subsref(obj,S);
                ret{e} = class(ele);                
            end
        end
    end
    
    
    methods (Static)
        %%%%%%%%%%%%%%%%%
        % to point set
        function [pS] = toPointSet(pL)
            pS = {};
            for e = 1:size(pL.d,2)
                pS{e} = phytoPoint(pL.d(:,e));
            end
        end
        
        
    end
        
end

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cat overload
[ret] = vertcat(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cat overload
[ret] = horzcat(varargin);
%}

