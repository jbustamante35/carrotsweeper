classdef typeConstraint < constraint

    properties
        allowList;
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = typeConstraint()
            obj.allowList = {};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add elements to allowable list
        function [] = addAllow(obj,type)
            if isa(type,'char')
                obj.allowList{end+1} = type;
            elseif isa(type,'cell')
                for e = 1:numel(type)
                    obj.allowList{end+1} = type{e};
                end
            end 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add elements to allowable list     
        function [b] = isAllow(obj,ele)
            b = typeConstraint.myisa(ele,obj.allowList);
        end
        
        
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % myisa
        function [r] = myisa(ele,type)
            if iscell(type)
                r = typeConstraint.wildCard(type);
                for e = 1:numel(type)
                    r(e+1) = builtin('isa',ele,type{e});    
                end
                r = any(r);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wild card type
        function [r] = wildCard(type)
            r = any(strcmp(type,'*'));
        end
    end
end
