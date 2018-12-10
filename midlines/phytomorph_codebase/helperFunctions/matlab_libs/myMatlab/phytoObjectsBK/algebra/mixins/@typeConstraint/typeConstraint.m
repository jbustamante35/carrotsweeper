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
            b = isa(ele,obj.allowList);
        end
        
    end
end
