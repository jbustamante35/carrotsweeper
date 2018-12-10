classdef mListable < matlab.mixin.Copyable

    properties
        mDi;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = mListable()
            obj.mDi = dictionary();
            % record methods
            met = methods(obj);
            for e = 1:numel(met)
                obj.mDi.insertTerm(met{e});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % isMethod
        function [ret] = isMethod(obj,mName)
            ret = isTerm(obj.mDi,mName);
        end
    end
    
    
end
