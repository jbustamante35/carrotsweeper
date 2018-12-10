classdef clusterForest < handle
    properties
        uniqueID;
        
        johnnyAppleseed;

        treeSet;

        extractFunc;
        clusterFunctionFunctionGenerator;
    end

    methods

        function [obj] = clusterForest(johnnyAppleseed)
            % assign unique ID
            obj.uniqueID = char(java.util.UUID.randomUUID);
            
            % assign pointer to arborist
            obj.johnnyAppleseed = johnnyAppleseed;
            % attach the arborist to this cluster forest
            obj.johnnyAppleseed.attachClusterForest(obj);

            % assign the extract function from john
            obj.extractFunc = johnnyAppleseed.extractFunc;

            % assign the cluster function fucntion generator from John
            obj.clusterFunctionFunctionGenerator = johnnyAppleseed.clusterFunctionFunctionGenerator;

            % init tree set
            obj.treeSet = {};
        end

        function [] = addTree(obj,tree)
            obj.treeSet{end+1} = tree;
        end

        function [tree] = getTree(obj,index)
            tree = obj.treeSet{index};
        end

        function [n] = numberOfTrees(obj)
            n = numel(obj.treeSet);
        end

        function [k p] = clusterImage(obj,data)
            for tree = 1:numel(obj.treeSet)
                if nargout == 1
                    [k(:,:,tree)] = obj.treeSet{tree}.clusterImage(data);
                elseif nargout == 2
                    [k(:,:,tree),p{tree}] = obj.treeSet{tree}.clusterImage(data);
                end
            end
        end

        function [k p] = clusterData(obj,data)
            for tree = 1:numel(obj.treeSet)
                if nargout == 1
                    [k(:,tree)] = obj.treeSet{tree}.cluster(data);
                elseif nargout == 2
                    [k(:,tree),p{tree}] = obj.treeSet{tree}.cluster(data);
                end
            end
        end
        
        function [] = pruneStumps(obj)
            for e = 1:obj.numberOfTrees
                if obj.treeSet{e}.numberOfLeaves == 1
                    rm(e) = true;
                else
                    rm(e) = false;
                end
            end
            obj.treeSet(rm) = [];
        end
   





        

        
    end
end