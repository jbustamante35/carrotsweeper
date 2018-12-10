classdef my_hc < myiF
    
    % my cluster class
    properties
        
    end
    
    % methods
    methods
        
        % constructor
        function [obj] = my_hc(varargin)            
            obj@myiF(@(x,y)my_hc.f_s(x,y),@(x)my_hc.i_s(x));
            if nargin == 1
                obj.p = varargin{1};
            end
        end
        
        % add cluster method        
        function [obj] = addMethod(obj,f,idx)
            if nargin == 2
                idx = numel(p) + 1;
            end
            p{idx} = f;
        end

    end
    
    methods(Static)
        
        % parameters
        function [ret] = i_s(d)
            fprintf(['Nope.not yet invertable\n']);
        end
        
        % parameters
        function [ret] = f_(d,p)
            sz = size(d.d,1);
            ret = ones(sz,1);            
            ret = hLabel(d.d,ret,p);            
        end
    end
    
end