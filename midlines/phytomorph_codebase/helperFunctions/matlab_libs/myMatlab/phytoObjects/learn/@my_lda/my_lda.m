classdef my_lda < myiF
    
    % my tensor class
    properties
        
    end
    
    % methods
    methods
        % constructor
        function [obj] = my_lda()
            obj@myiF(@(d,p)my_lda.f_s(d,p),@(d,c,p)my_lda.i_s(d,c,p));
        end
    end
    
    methods(Static)
        
        % parameters
        function [ret] = i_s(d,c,p)
            ret = myLDA(d,c);
        end
        
        % parameters
        function [ret] = f_s(d,p)
            ret = d.domain*p;
        end
    end
    
end