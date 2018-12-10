classdef clusterNode < handle
    
    
    properties
        uniqueID;
        
        featureIDX;
        idxSelectorFunction;

        level;
        selfIndex;
        parentIndex;
        treePointer;
        maxClusterLevels;
        childList;
        
        clusterFunction;
        levelClusterFunction;
        probFunction;


        % generates the cluster parameters to send to cluster function
        clusterFunctionFunctionGenerator;
        suggestNumberOfClustersFunction;
    end
    
    
    methods
        function [obj] = clusterNode(treePointer,parentIndex,level,suggestNumberOfClustersFunction,clusterFunctionFunctionGenerator,idxSelectorFunction)
            obj.uniqueID = char(java.util.UUID.randomUUID);
            obj.childList = {};
            obj.treePointer = treePointer;
            obj.parentIndex = parentIndex;
            obj.level = level;
            obj.suggestNumberOfClustersFunction = suggestNumberOfClustersFunction;
            obj.clusterFunctionFunctionGenerator = clusterFunctionFunctionGenerator;
            obj.idxSelectorFunction = idxSelectorFunction;
        end
        
        function [] = build(obj,data)
            try
                
                % suggest feature(s)
                obj.featureIDX = obj.idxSelectorFunction(data,obj.level);



                % suggest number of clusters
                obj.maxClusterLevels = obj.suggestNumberOfClustersFunction(data(:,obj.featureIDX), obj.level);
                obj.maxClusterLevels = obj.maxClusterLevels + 1;
                UQ = [];
                
                while numel(UQ) ~= obj.maxClusterLevels

                    
                    % if not all clusters are used
                    if numel(UQ) ~= obj.maxClusterLevels
                        obj.maxClusterLevels = obj.maxClusterLevels - 1;
                    end
                    
                    % if one cluster than always cluster with 1 else build cluster func
                    if obj.maxClusterLevels == 1
                        func = @(data)obj.selfIndex;
                        % set cluster function
                        obj.setClusterFunction(func);
                        % report self as leaf
                        obj.treePointer.insertLeaf(obj);
                        % set prob to one
                        obj.setProbFunction(@(X)ones(size(X,1),1));
                    else

                        gmm = obj.clusterFunctionFunctionGenerator(data(:,obj.featureIDX),obj.maxClusterLevels);
                        % create cluster function
                        func = @(X)cluster(gmm,X);
                        % set cluster function
                        obj.setClusterFunction(func);
                        % set prob function
                        obj.setProbFunction(@(X)posterior(gmm,X));

                    end

                    % set the level cluster function to return the unique ID of node
                    obj.levelClusterFunction = @(data)obj.selfIndex;

                    % cluster the data
                    k = obj.cluster(data);
                    
                    % get number of populated clusters
                    UQ = unique(k);
                    
                    
                end

                % attempt to subcluster if obj.maxClusterLevels ~= 1
                if obj.maxClusterLevels ~= 1
                    for level = 1:obj.maxClusterLevels
                        sidx = find(k == level);
                        if isempty(sidx)
                            stop = 1;
                        end
                        childNode = clusterNode(obj.treePointer,obj.selfIndex,obj.level+1,obj.suggestNumberOfClustersFunction,obj.clusterFunctionFunctionGenerator,obj.idxSelectorFunction);
                        % must insert before build 
                        obj.treePointer.insertNode(childNode);
                        childNode.build(data(sidx,:));
                        obj.attachChild(childNode);
                    end
                end
            catch ME
                ME
            end
            
        end
        
        function [] = setClusterFunction(obj,clusterFunction)
            obj.clusterFunction = clusterFunction;
        end

        function [] = setProbFunction(obj,probFunction)
            obj.probFunction = probFunction;
        end
    
        function [p] = prob(obj,data)
            p = obj.probFunction(data);
        end

        function [k,pi] = cluster(obj,data,k,idx,pi,pstate)

            if nargin == 2
                k = zeros(size(data,1),1);
                idx = (1:size(data,1))';
            end

            ki = obj.clusterFunction(data(:,obj.featureIDX));

            if nargout == 2
                pnew = obj.probFunction(data(:,obj.featureIDX));
            end

            % if there are children
            if ~isempty(obj.childList)
                
                if nargout == 2
                    
                    
                    for e = 1:numel(obj.childList)
                        %fidx = find(ki == e);
                        %fidx = 1:numel(ki);
                        pupdate = pnew(:,e).*pstate;
                        [k,pi] = obj.childList{e}.cluster(data,k,idx,pi,pupdate);
                    end
                    
                    
                else
                    for e = 1:numel(obj.childList)
                        fidx = find(ki == e);
                        [k] = obj.childList{e}.cluster(data(fidx,:),k,idx(fidx));
                    end
                end

            else
                
                
                if nargout == 1
                    k(idx) = ki;
                elseif nargout == 2
                    index = obj.treePointer.findLeafNodeIndex(obj);
                    pi(idx,index) = pnew.*pstate;
                end
            end


        end
        
        function [] = attachChild(obj,childNode)
            obj.childList{end+1} = childNode;
        end
        
    end
end