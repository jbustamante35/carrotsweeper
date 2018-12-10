classdef myTree < handle
    properties
        nodelist = {};
        nodeIDlist = [];
        numNodes;
        headnode;
    end
    
    methods
        function [obj] = myTree()
            obj.numNodes = 0;
            obj.headnode = tree_node();
            obj.registerNodeWithTree(obj.headnode);
        end
        
        function [] = attachChildToParent(obj,parent,child)
            obj.registerNodeWithTree(child);
            % attach child to parent
            parent.attachChildren(child);
        end
        
        function [hn] = getHeadNode(obj)
            hn = obj.headnode;
        end
        
        function [curves] = findSpecial(obj)            
            idx = [];
            cnt = 1;
            for e = 1:numel(obj.nodelist)
                nc = obj.nodelist{e}.getNumberOfChildren();
                ns = obj.nodelist{e}.getNumberOfSibs();
                if ((ns == 0 & nc > 1) | (ns == 0 & nc == 0)) & e ~= 1
                    curves(cnt) = obj.nodelist{e}.parentlist{1}.data;
                    cnt = cnt + 1;
                end
            end
            if cnt == 1;
                curves = [];
            end
        end
        
        
    end
    
    methods (Access = private)
        function [] = registerNodeWithTree(obj,node)
            obj.numNodes = obj.numNodes + 1;
            obj.nodeIDlist(obj.numNodes) = node.id;
            obj.nodelist{obj.numNodes} = node;
        end        
    end
end