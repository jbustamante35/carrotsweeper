classdef myHS_2 < myHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % abstract def of hetero-geneous container
    % two for two types - note that here that one
    % of the types is the container
    properties       
        eType = '';         % element type
        extractFunction;    % extract function
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = myHS_2(varargin)
            % is non-default to nothing useful
            if nargin == 0
                obj.S = {};
                obj.eType = '';
                obj.extractFunction = '';
            end
            % if one then should be element type
            if nargin >= 1
                % set the element type
                obj.eType = varargin{1};
                % init the set
                obj.S = {};
                % type operator for extraction
                obj.extractFunction = typeOp(obj.eType,class(obj));
            end
            %
            if nargin == 2
                obj.S = varargin{1};
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [b] = allow(obj,e)
            
            % if isa allowed element type
            b = isa(e,obj.eType);
            
            % is proper container
            if isa(e,class(obj))
                % search container's elements
                for i = 1:numel(e)
                    % build of subs ref for internal
                    subs.type = '{}';
                    subs.subs{1} = i;
                    % call to subs ref
                    c = subsref(e,subs);
                    % ask if each element is allowed
                    b(i) = obj.allow(c);
                end
                b = all(b==1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % putElement
        function [] = putElement(obj,e,i)
            if nargin == 2
                i = numel(obj.S)+1;
            end
            S(1).type = '{}';
            S(1).subs{1} = i;
            obj.subsasgn(S,e);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getElement
        function [r] = getElement(obj,i)            
            S(1).type = '{}';
            S(1).subs{1} = i;
            r = obj.subsref(S);
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [ret] = subsasgn(obj,S,B)
            ret = obj;
            if obj.allow(B);
                obj = subsasgn@myHS(obj,S,B);
            else
                fprintf(['object@class: ' class(B) ' was rejected\n']);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % harvest
        function [r] = harvest(obj,op)
            % eval harvest over graph
            op.op(obj);
            % return results
            r = op.r;
            % clear res
            op.clearH();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % default view implementation
        function [] = view(obj,h,varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % persist any props through function call            
            obj.view_props.props = viewable.parseView(varargin,obj.view_props.props);        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % harvest first level nodes
            obj.extractFunction.setN(1); 
            T = obj.harvest(obj.extractFunction);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % view each node - each node (default myT, must have view
            % defined)
            for e = 1:numel(T)
                t = T{e};
                t = t.getElement(1);
                t.setView(obj.view_props);
                t.view(h);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obtain representation
        function [R] = rep(obj)
            R = [];    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % harvest first level nodes
            obj.extractFunction.setN(1); 
            T = obj.harvest(obj.extractFunction);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % view each node - each node (default myT, must have view defined)
            for e = 1:numel(T)
                t = T{e}.subset(1);
                R = [R t.d];
            end
            R = mean(R,2);
        end
        
    end
end
%{
    S = myHS_2('char');
    S{1} = 'hello';
    S{2} = 'world';
%}
