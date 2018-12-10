classdef my_lsf < myiF
    %%%%%%%%%%%%%
    % "invert" functions
    properties
    end
    
    %%%%%%%%%%%%%
    % methods
    methods
        % constructor
        function [obj] = my_lsf()
            obj@myiF(@(d,p)my_lsf.f_s(d,p),@(d,c,p)my_lsf.i_s(d,c,p));
        end
    end

    methods(Static)

        % parameters
        function [ret] = i_s(d,c,p)
            ret = mldivide(d.d,c.d);            
        end
        
        % parameters
        function [ret] = f_s(d,p)
            ret = d.d*p;
            ret = myT(ret);
        end
    end
    
end