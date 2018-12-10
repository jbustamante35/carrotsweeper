classdef johnMuir < handle
    
    properties
        metricFunc
        overlayFunc
        oPath = '';
    end
    
    methods
        function [obj] = johnMuir()
            metricFunc = [];
            overlayFunc = [];
        end
        
        function [] = attachMetric(obj,metricFunc)
            obj.metricFunc = metricFunc;
        end
        
        function [] = attachOverlayFunction(obj,overlayFunc,oPath)
            obj.overlayFunc = overlayFunc;
            obj.oPath = oPath;
        end
        
        function [sT] = measureTree(obj,cT,fileList,maxSubGroups)
            try
                uniqueLabels = cell2mat(cT.leafNodeIndexList);
                subReClustering = maxSubGroups;
                sT = {};
                for f = 1:numel(fileList)
                    iT = table;
                    % get orginal image
                    oI = double(imread(fileList{f}));
                    % get file parts
                    [p,nm,ext] = fileparts(fileList{f});
                    % cluster image
                    k = cT.clusterImage(fileList{f});

                    
                    tmpP = [obj.oPath filesep cT.uniqueID filesep];
                    mkdir(tmpP)

                    for c = 1:subReClustering
                        groupings = nchoosek(uniqueLabels,c);
                        for g = 1:size(groupings,1)

                            columnHeader = ['c' num2str(c) '_g' num2str(g)];
                            
                            
                            % generate mask
                            mask = zeros(size(k));
                            for l = 1:size(groupings,2)
                                mask = mask | k == groupings(g,l);
                            end
                            % measure mask
                            m = obj.metricFunc(mask);
                            
                            
                            % START HERE _ NEED UNIQUE NAME FOR TREE
                            % overlay is available
                            if ~isempty(obj.overlayFunc)
                                
                                oName = [tmpP ' grp_' columnHeader filesep];
                                mkdir(oName);
                                oName = [oName '{treeID_' cT.uniqueID '}{fileName_' nm '}{grp_' columnHeader '}.tif'];
                                obj.overlayFunc(oI,mask,oName);
                            end
                            
                            
                            
                            rowNames = m.Properties.VariableNames;
                            for r = 1:numel(rowNames)
                                iT{rowNames{r},columnHeader} = m{1,rowNames{r}};
                            end
                        end
                    end
                    sT{f} = iT;
                end
            catch ME
                ME
            end
            
        end

        function [m] = measureForest(obj,forest,fileList,maxSubGroups,treeList)
            if nargin ==4
                treeList = 1:forest.numberOfTrees();
            end
            for e = treeList
                m{e} = obj.measureTree(forest.getTree(e),fileList,maxSubGroups);
            end
        end
    end
end