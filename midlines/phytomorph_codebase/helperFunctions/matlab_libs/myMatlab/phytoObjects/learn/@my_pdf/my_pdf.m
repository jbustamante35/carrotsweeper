classdef my_pdf < myiF
    
    % my tensor class
    properties
        
    end
    
    % methods
    methods
        % constructor
        function [obj] = my_pdf()
            obj@myiF(@(d,p)my_pdf.f_s(d,p),@(d,c,p)my_pdf.i_s(d,c,p));
        end
    end
    
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%
        % parameters
        function [ret] = i_s(d,c,p)
            UQ = unique(c);
            for u = 1:numel(UQ)
                idx = find(c == UQ(u));
                mu = mean(d.domain(idx,:),1);
                v = var(d.domain(idx,:),1);
                ret = [ret;[mu;v]];
            end
        end
        
        %%%%%%%%%%%%%%%%%%
        % parameters
        function [ret] = f_(d,p)
            sz = size(p,2);
            for f = 1:size(p,1)            
                ret(f) = mvnpdf(d,p(f,1:sz/2),p(f,sz/2+1:end));
            end
        end
    end
    
end