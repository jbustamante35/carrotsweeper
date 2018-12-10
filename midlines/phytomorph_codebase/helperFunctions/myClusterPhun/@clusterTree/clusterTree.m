classdef clusterTree < handle
    properties
        
        uniqueID;
        
        extractFunc;
        
        idxSelectorFunction;

        suggestNumberofClusters;
        generateClusterFunction;
        clusterNodeList;
        
        leafNodeIndexList;
        
        suggestNumberOfClusterFunction;
        clusterFunctionFunctionGenerator;
        
        
    end
    
    methods
        function [obj] = clusterTree(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc,idxSelectorFunction)
            obj.uniqueID = char(java.util.UUID.randomUUID);
            
            obj.idxSelectorFunction = idxSelectorFunction;

            % this function suggests the number of branches. 
            % must take the level and data into account.
            % this can be created by an arborist.
            obj.suggestNumberOfClusterFunction = suggestNumberOfClusterFunction;
            

            obj.clusterFunctionFunctionGenerator = clusterFunctionFunctionGenerator;
            obj.extractFunc = extractFunc;
            obj.clusterNodeList = {};
            obj.leafNodeIndexList = {};
        end
        
        function [] = build(obj,data)
            headNode = clusterNode(obj,0,1,obj.suggestNumberOfClusterFunction,obj.clusterFunctionFunctionGenerator,obj.idxSelectorFunction);
            obj.insertNode(headNode);
            headNode.build(data);
        end
        
        function [] = insertNode(obj,node)
            obj.clusterNodeList{end+1} = node;
            node.selfIndex = numel(obj.clusterNodeList);
        end
        
        function [k,p] = cluster(obj,data)

            idx = (1:size(data,1))';
            k = zeros(size(data,1),1);
            p = zeros(size(data,1),obj.numberOfLeaves());
            pstate = ones(size(data,1),1);
            if nargout == 1
                k = obj.clusterNodeList{1}.cluster(data,k,idx);
            elseif nargout == 2
                [k,p] = obj.clusterNodeList{1}.cluster(data,k,idx,p,pstate);
            end
        end
        
        function [] = insertLeaf(obj,leafNode)
            obj.leafNodeIndexList{end+1} = leafNode.selfIndex;
        end

        function [index] = findLeafNodeIndex(obj,leafNode)
            index = [];
            try
                if isa(leafNode,'clusterNode')
                    searchIndex = leafNode.selfIndex;
                else
                    searchIndex = leafNode;
                end
                for l = 1:numel(obj.leafNodeIndexList)
                    if (obj.leafNodeIndexList{l} == searchIndex)
                        index = l;
                    end
                end
                if isempty(index)
                    stop = 1;
                end
                
            catch ME
                ME;
            end
        end
        
        function [n] = numberOfLeaves(obj)
            n = numel(obj.leafNodeIndexList);
        end
        
        function [k,p] = clusterImage(obj,data)
            if ischar(data)
                data = double(imread(data));
            end
            
            sz = size(data);
            
            data = obj.extractFunc(data);
            
            
            if nargout == 1
                [k] = obj.cluster(data);
                k = reshape(k,sz(1:2));
                
                
                UQ = unique(k);
                newk = k;
                for u = 1:numel(UQ)
                    newLabel = obj.findLeafNodeIndex(UQ(u));
                    fidx = find(k == UQ(u));
                    newk(fidx) = newLabel;
                end
                k = newk;
                
            elseif nargout == 2
                [k,p] = obj.cluster(data);
                k = reshape(k,sz(1:2));
                p = reshape(p,[sz(1:2) obj.numberOfLeaves]);
                [~,k] = max(p,[],3);
            end
            
            
            
        end
    end
end

%{


%}