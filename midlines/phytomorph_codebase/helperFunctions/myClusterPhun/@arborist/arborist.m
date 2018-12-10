classdef arborist < handle
    properties
        % parameters for generateing planting tree function
        maxDepth;
        maxBranchNumber;
        % function to call and leave open for data and level
        % will suggest the number of branches on the tree based
        % on data and levels
        functionGenerator;

        % IDX selector function
        idxSelectorFunction;

        % clusterForest to plant - current forest
        cF;

        % generates the cluster parameters to send to cluster function
        clusterFunctionFunctionGenerator;


        % feature extraction function
        extractFunc;
    end

    methods
        function [obj] = arborist(functionGenerator,clusterFunctionFunctionGenerator,extractFunc,maxDepth,maxBranchNumber,idxSelectorFunction)
            % for generating tree growth parameters
            obj.functionGenerator = functionGenerator;
            obj.maxDepth = maxDepth;
            obj.maxBranchNumber = maxBranchNumber;

            % for planting trees
            obj.clusterFunctionFunctionGenerator = clusterFunctionFunctionGenerator;
            obj.extractFunc = extractFunc;

            obj.idxSelectorFunction = idxSelectorFunction;
        end

        function [] = attachClusterForest(obj,cF)
            obj.cF = cF;
        end

        function [func] = suggestATreeSuggestNumberOfClusterFunction(obj)
            
            treePara(1) = randi(obj.maxDepth);
            treePara(2) = randi(obj.maxBranchNumber);

            func = @(data,level)obj.functionGenerator(data,level,treePara);
        end


        function [forest] = plantTrees(obj,data,numberOfTrees)
            forest = clusterForest(obj);
            for t = 1:numberOfTrees
                treeBuildingFunction = obj.suggestATreeSuggestNumberOfClusterFunction;
                tree = clusterTree(treeBuildingFunction,obj.clusterFunctionFunctionGenerator,obj.extractFunc,obj.idxSelectorFunction);
                tree.build(data);
                forest.addTree(tree);  
            end
            forest.pruneStumps();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data :- data to sample
        % sampleParameters{1} :- number of sets to sample
        % sampleParameters{2} :- number of objects per set to sample
        % sampleParameters{3} :- number of elements per object
        % sampleParameters{4} :- resize parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [S] = sampleElements(obj,data,sampleParameters)
            
            if iscell(data)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % N : = number of sets to sample
                N = sampleParameters(1);
                % sN : = number of objects per set to sample
                sN = sampleParameters(2);
                % number of elements per object
                M = sampleParameters(3);
                % resize parameters
                re = sampleParameters(4);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if data is a charater - then read
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fileName = data{1};
                if ~ischar(fileName)
                    fileName = fileName{1};
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % run extract function 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tmp = obj.extractFunc(imread(fileName));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % pre-allocate the sample data
                %S = zeros(N*sN*M,size(tmp,2));
                toSample = randi(numel(data),1,N);
                str = 1;
                stp = str + M - 1;
                
                cnt = 1;
                % for 1:N
                for s = 1:N
                    % 
                    tmpToSample = data{toSample(s)};
                    if ischar(tmpToSample)
                        tmpToSample = {tmpToSample};
                    end
                    subS = randi(numel(tmpToSample),1,sN);
                    % over the sample number in the set
                    for ss = 1:numel(subS)
                        randomList{cnt} = tmpToSample{subS(ss)};
                        cnt = cnt + 1;
                    end
                   
                end
                
                
                
                % over the sample number in the set
                parfor ss = 1:numel(randomList)
                    fprintf(['start sampling:' num2str(s) ':' num2str(ss) ':' num2str(N) '\n']);
                    % read the image
                    I = imread(randomList{ss});
                    % resize if needed
                    if re ~= 1
                        I = imresize(I,re);
                    end
                    % run extract function
                    I = obj.extractFunc(I);
                    % obtain random schuffle
                    idx = randperm(size(I,1));
                    % schuffle
                    I = I(idx,:);
                    % sample and store
                    S{ss} = I(idx(1:M),:);
                    %str = stp + 1;
                    %stp = str + M - 1;
                    fprintf(['stop sampling:' num2str(s) ':' num2str(ss) ':' num2str(N) '\n']);
                end
                S = cell2mat(S');
            else
                N = sampleParameters(1);
                M = sampleParameters(2);
                tmp = obj.extractFunc(data(:,:,:,1));
                S = zeros(M*N,size(tmp,2));
                
               

                toSample = randi(size(data,4),1,N);
                str = 1;
                stp = str + M - 1;
                for s = 1:N
                    fprintf(['start sampling:' num2str(s) ':' num2str(N) '\n']);
                    tmp = obj.extractFunc(data(:,:,:,toSample(s)));
                    idx = randperm(size(tmp,1));
                    tmp = tmp(idx,:);
                    S(str:stp,:) = tmp(idx(1:M),:);
                    str = stp + 1;
                    stp = str + M - 1;
                    fprintf(['end sampling:' num2str(s) ':' num2str(N) '\n']);
                end
            end
            
            
            
            
        end


    end
end