classdef tree_node < handle

    properties
        id;
        childlist = {};
        parentlist = {};
        data;
    end
    
    methods
        function [obj] = tree_node()
            obj.id = string2hash([datestr(now) num2str(rand(1,100))]);
        end 
        
        function [] = attachChildren(obj,child)
            obj.childlist{end+1} = child;
            child.parentlist{end+1} = obj;
        end
        
        function [] = attachData(obj,data)
            obj.data = data;
        end
        
        function [treeWidth] = getTreeWidth(obj)            
            curCount = obj.getChildrenNumber([],1);
            UQ_depth = unique(curCount(:,1));
            for u = 1:numel(UQ_depth)
                fidx = find(curCount(:,1)==UQ_depth(u));
                treeWidth(u) = sum(curCount(fidx,2));
            end
        end
        
        function [n] = getNumberOfChildren(obj)
            n = numel(obj.childlist);
        end
        
        function [n] = getNumberOfSibs(obj)
            if ~isempty(obj.parentlist)
                n = numel(obj.parentlist{1}.childlist)-1;
            else
                n = 0;
            end
        end
        
    end
    
    methods (Access = private)
        function [curCount] = getChildrenNumber(obj,curCount,curDepth)
            curCount = [curCount;[curDepth+1 numel(obj.childlist)]];
            for e = 1:numel(obj.childlist)
                curCount = obj.childlist{e}.getChildrenNumber(curCount,curDepth+1);
            end
        end
    end
    
    
    
end