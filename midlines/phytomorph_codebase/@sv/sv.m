classdef sv < handle
    properties
        v;
    end
    
    methods
        
        function [obj] = sv(sz1,sz2)
            obj.v = sparse(sz1,sz2);
        end
        
        function [w] = mtimes(obj,v)
            w = conj(obj.v)'*conj(v.v);
        end
        
        function [B] = subsref(obj,S)
            if strcmp(S.type,'.')
                B = obj.v;
            else
                B = subsref(obj.v,S);
            end
            
        end
        
        function [obj] = subsasgn(obj,S,A)
            obj.v = subsasgn(obj.v,S,A);
        end
        
        function [sz] = size(obj,varargin)
            sz = size(obj.v,varargin{:});
        end
        
        function [b] = issparse(obj)
            b = issparse(obj.v);
        end
    end
end