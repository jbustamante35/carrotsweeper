classdef my_cm < myiF
    
    % my cluster mapper
    properties
        
    end
    
    % methods
    methods
        
        % constructor
        function [obj] = my_cm(varargin)            
            obj@myiF(@(d,p)my_cm.f_s(d,p),@(d,c,p)my_cm.i_s(d,c,p));
            if nargin == 1
                obj.ip = varargin{1};
            end
        end
        
        % add cluster method        
        function [obj] = addMethod(obj,f,idx)
            if nargin == 2
                idx = numel(obj.ip) + 1;
            end
            obj.ip{idx} = f;
        end

    end
    
    methods(Static)
        
        % parameters
        function [cc] = i_s(d,c,p)
            % cluster
            sz = size(d.d,1);
            ret = ones(sz,1);            
            ret = hLabel(d.d,ret,p);
            % obtain means of cluster centers
            UQ = unique(ret);
            for u = 1:numel(UQ)
                idx = find(ret == UQ(u));
                cc(u,:) = mean(d.d(idx,:));
            end
        end
        
        % parameters
        function [ret] = f_s(d,p)
            % for each data point
            for e = 1:size(d.d,1)
                delta = bsxfun(@minus,p,d.d(e,:));
                delta = sum(delta.*delta,2);
                [JUNK,ret(e)] = min(delta);
            end       
        end
    end
    
end