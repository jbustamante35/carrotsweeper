classdef dictionary < myHS_X
    
    properties
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = dictionary(varargin)
            obj@myHS_X('char');
            if nargin == 1
                for e = 1:numel(varargin{1})
                    % assign
                    S(1).type = '{}';
                    S(1).subs{1} = e;
                    subsasgn(obj,S,varargin{1}{e});
                end
            end
        end
        %%%%%%%%%%%%%%%%%%
        % insert
        function [] = insertTerm(obj,term,n)
            if nargin == 2
                n = numel(obj.S)+1;
            end
            obj.S{n} = term;
        end
        %%%%%%%%%%%%%%%%%%
        % lookupTerm
        function [key] = lookUp(obj,term)
            key = find(strcmp(obj.S,term));
        end
        %%%%%%%%%%%%%%%%%%
        % isTerm
        function [ret] = isTerm(obj,key)
            ret = ~isempty(lookUp(obj,key));
        end
        %%%%%%%%%%%%%%%%%%
        % matchTerm
        function [ret] = seachTerm(obj,sKey)
            ret = zeros(1,numel(obj));
            for e = 1:numel(obj)
                sT = obj.S{e};
                ret(e) = strcmp(sT(1:numel(sKey)),sKey);
            end
            ret = find(ret);
        end
        
        
    end
end