classdef myHS_X < myHS & typeConstraint
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % abstract def of hetero-geneous container
    % one for one type
    properties       
        extractFunction;    % extract function
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = myHS_X(varargin)
            % call super constructor
            obj = obj@myHS();
            obj = obj@typeConstraint();
            
            % is non-default to nothing useful
            if nargin == 0
                obj.extractFunction = '';
            end
            
            % if one then should be element type
            if nargin >= 1
                % set the element type
                obj.addAllow(varargin{1});
                % type operator for extraction
                %obj.extractFunction = typeOp(obj.eType,class(obj));
            end
            
            % init the elements of the set S
            if nargin == 2
                fprintf(['problem here.\n']);
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
            if obj.isAllow(B);
                ret = subsasgn@myHS(obj,S(1),B);
                S(1) = [];
            else
                fprintf(['object@class: ' class(B) ' was rejected\n']);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [r] = subsref(obj,S)
            try
                if nargout == 0
                   subsref@myHS(obj,S);
                else
                    r = subsref@myHS(obj,S);
                end
            catch ME
                try
                    F = str2func(['@(varargin)' S(1).subs '(varargin{:})']);
                    G = {F S(2).subs{:}};                    
                    % pass along other args
                    if nargout == 0
                        obj.distrib(G{:});
                    else
                        r = obj.distrib(G{:});
                        %r = obj.distrib(F,S);
                    end 
                catch ME2
                    ME2
                end
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
        % obtain representation
        function [] = rep(obj)
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view
        function [] = view(obj,h,vProps)
            S(1).type = '{}';        
            for e = 1:numel(obj)
                S(1).subs{1} = e;
                ele = subsref(obj,S);
                ele.view(h,vProps);
                if isa(ele,'phytoApoint')
                    d = ele(1);
                    text(d(1),d(2),num2str(e));
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % type of elements in set - first level deep
        function [r] = distrib(obj,varargin)
            op = varargin{1};
            varargin(1) = [];
            if nargout == 1;r = myHS_X('*');end
            if isa(op,'char');op=str2func(op);end
            S(1).type = '{}';
            for e = 1:numel(obj)
                S(1).subs{1} = e;
                if nargout == 1
                    ele = subsref(obj,S);
                    G = {ele varargin{:}};
                    ele = op(G{:});
                    r.putElement(ele);
                else
                    ele = subsref(obj,S);
                    G = {ele varargin{:}};
                    op(G{:});
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % type of elements in set - first level deep
        function [r] = par_distrib(obj,op)
            if nargout == 1;r = myHS_X('*');end
            if isa(op,'char');op=str2func(op);end
            S(1).type = '{}';
            for e = 1:numel(obj)
                S(1).subs{1} = e;
                if nargout == 1
                    ele = subsref(obj,S);
                    ele = op(ele);
                    r.putElement(ele);
                else
                    op(subsref(obj,S));
                end
            end
        end
    end
end



