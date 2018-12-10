classdef myHS_1 < myHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % abstract def of hetero-geneous container
    % one for one type
    properties       
        eType = '';         % element type
        extractFunction;    % extract function
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = myHS_1(varargin)
            
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
            
            % init the elements of the set S
            if nargin == 2
                obj.S = varargin{2};
            end
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [b] = allow(obj,e)
            % if isa allowed element type
            b = isa(e,obj.eType);
        end
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [ret] = subsasgn(obj,S,B)
            ret = obj;
            if obj.allow(B);
                ret = subsasgn@myHS(obj,S,B);
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
        % obtain representation
        function [] = rep(obj)
            
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view
        function [] = view(obj,h,varargin)            
            S(1).type = '{}';        
            for e = 1:numel(obj)
                S(1).subs{1} = e;
                ele = subsref(obj,S);
                ele.view(h);
            end
        end
        
    end
end



