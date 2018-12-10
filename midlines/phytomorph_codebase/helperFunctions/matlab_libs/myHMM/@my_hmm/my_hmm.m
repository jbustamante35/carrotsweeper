classdef my_hmm < handle

    properties
        state = [];
        NodeList = {};
        dn = [];
        %history = 20;
        %history = 5;
        history = 150;
        %history = 100;
        %history = 400;
        %history = 1000;
        DYN = 1;
        needPRE = 1;
        TStore = [];
    end
    
    
    methods
        
        function [obj] = my_hmm()
        
        end
        
        function [] = addNode(obj,node)
            obj.NodeList{end+1} = node;
        end
        
        function [emm] = getObservations(obj,N)
            nn = obj.NodeList{obj.state}.getNextNode();
            for e = 1:N-1          
                nn = nn.getNextNode();
                emm(:,e) = nn.observe();
            end            
        end
        
        function [prob] = forwardCalc(obj,observation,prob)
            for o = 1:size(observation,2)
                for e = 1:numel(obj.NodeList)
                    tmpProb = obj.NodeList{e}.forwardCalc(obj.NodeList,observation(:,o),obj.dn);
                    P(e) = sum(prob.*tmpProb);
                end
                prob = P;
            end
        end
        
        %{
        function [tm] = buildTransitionMatrix(obj)
            for source = 1:numel(obj.NodeList)
                for target = 1:numel(obj.NodeList)
                    tm(source,target) = obj.NodeList{source}.toTransitionProb(obj.NodeList{target});
                end                
            end
        end
        %}

        function [tm] = buildTransitionMatrix(obj,stateObservation)
            tm = [];
            for source = 1:numel(obj.NodeList)
                sourceNode = obj.NodeList{source};
                sourceNodeNumber = obj.translateNodeToNumber(sourceNode);
                for target = 1:numel(obj.NodeList)
                    targetNode = obj.NodeList{target};
                    tm(source,target) = sourceNode.getTransitionProbToNode(targetNode,stateObservation,{sourceNodeNumber});
                end
            end
        end
        
        function [states,prob] = Viterbi(obj,observation,observationLabel,forceStartState)
            
            
            
            
            states = [];
            
            % init paths
            for e = 1:numel(obj.NodeList)
                paths{e} = [];
            end
            
            % init cost per path
            for n = 1:numel(obj.NodeList)
                costPerPath(1,n) = log(obj.NodeList{n}.getObservationProb(observation(:,1),observationLabel,obj.dn));
            end
            
            if nargin == 4
                costPerPath = -inf*ones(size(costPerPath));
                costPerPath(forceStartState) = 0;
            end
            
            historyCost = [];
            for o = 2:size(observation,2)
                % select the path
                [paths,costPerPath] = obj.Viterbi_choice([],costPerPath,paths);
                % update the node cost by the emission amount
                for n = 1:numel(obj.NodeList)
                    emissionProb(n) = log(obj.NodeList{n}.getObservationProb(observation(:,o),observationLabel,obj.dn));
                    %emissionProb(n) = (obj.NodeList{n}.getObservationProb(observation(:,o),observationLabel,obj.dn));
                end
                %{
                % the below mixture may have been added to force a state -
                % it has been removed to fix the maize code
                if o == size(observation,2)
                    costPerPath = .01*costPerPath + .99*emissionProb;
                else
                    costPerPath = costPerPath + emissionProb;
                end
                %}
                PH{o} = paths;
                historyCost = [historyCost;costPerPath];
                costPerPath = costPerPath + emissionProb;
                %emissionProbLOG(o,:) = emissionProb;
                %costPerPathLOG(o,:) = costPerPath;
            end
            
            [prob idx] = max(costPerPath);
            states = paths{idx};
           
            
            
            
            
        end
        
        function [] = update(obj,updateData,observationLabels)
            
            
            transC = zeros(numel(obj.NodeList));
            stateC = zeros(1,numel(obj.NodeList));
            for e = 1:numel(updateData)
                [tmpT tmpS] = obj.countTransitions(updateData{e}(1,:),numel(obj.NodeList));
                transC = tmpT + transC;
                stateC = tmpS + stateC;
                fprintf(['Done counting @:' num2str(e) ':' num2str(numel(updateData)) '\n']);
            end
            stateC = stateC';
            stateC = sum(transC,2);
            
            newT = bsxfun(@times,transC,stateC.^-1);
            for source = 1:size(newT,1)
                for target = 1:size(newT,2)
                    if newT(source,target) ~= 0
                        obj.NodeList{source}.upDateTransitionProb(obj.NodeList{target},newT(source,target));
                    end
                end
            end
            
            
            
            
            
            UQ = unique(observationLabels);
            for s = 1:numel(obj.NodeList)
                tmpObservations = [];
                for e = 1:numel(updateData)
                    states = updateData{e}(1,:);
                    tmpObservations = [tmpObservations updateData{e}(2:end,states==s)];
                end
                for u = 1:numel(UQ)
                    sig = tmpObservations(observationLabels==UQ(u),:)';
                    %{
                    obj.NodeList{s}.Distribution{u}.mu = mean(sig,1);
                    obj.NodeList{s}.Distribution{u}.sigma = cov(sig);
                    %}
                    obj.NodeList{s}.Distribution{u}.update(sig);
                end
            end
                
            
        end
        
        function [idx] = detectIllegalTransitions(obj,data)
            filter = obj.buildTransitionMatrix();
            fidx = find(filter==0);
            idx = [];
            for e = 1:numel(data)
                [tmpT tmpS] = obj.countTransitions(data{e}(1,:),numel(obj.NodeList));
                if ~all(tmpT(fidx)==0)
                    idx = [idx e];
                end
            end
        end
        
        function [tmCounts stCounts] = countTransitions(obj,states,Tstates)
            % transistion counts
            %trans = [im2col(states,[1 2]) [states(end);states(1)]];
            trans = im2col(states,[1 2]);
            for i = 1:Tstates
                for j = 1:Tstates
                    sig = [i;j];
                    for k = 1:size(trans,2)
                        b(k) = all(trans(:,k)==sig);
                    end
                    tmCounts(i,j) = sum(b);
                end
            end
            % state counts
            for i = 1:Tstates
                stCounts(i) = sum(states==i);
            end
        end
    end
    
    methods (Access=private)
        function [paths,costPerPath] = Viterbi_choice(obj,transitionMatrix,costPerPath,paths)
            if obj.DYN
                % build up the trans matrix
                for target = 1:numel(obj.NodeList)
                    transCost = [];
                    targetNode = obj.NodeList{target};
                    for source = 1:numel(paths)
                        if ~isempty(paths{source})
                            sourceNumber = paths{source}(end);
                            len = numel(paths{source});
                            str = max(numel(paths{source})-obj.history,1);
                            historyToFeed = paths{source}(str:end);
                            %{
                            if all(historyToFeed == 4)
                                stop = 1;
                            end
                            %}
                            transCost(source) = log(obj.NodeList{sourceNumber}.getTransitionProbToNode(targetNode,historyToFeed,{sourceNumber len}));
                        else
                            len = numel(paths{source});
                            %transCost(source) = costPerPath(source); % question this
                            transCost(source) = log(obj.NodeList{source}.getTransitionProbToNode(targetNode,[],{target len}));
                        end
                    end
                    (transCost.*costPerPath.^-1);
                    [transitionCost(target) selPath(target)] = max(transCost+costPerPath);   
                end
            else
                if obj.needPRE
                    % build up the trans matrix
                    for target = 1:numel(obj.NodeList)
                        transCost = [];
                        targetNode = obj.NodeList{target};
                        for source = 1:numel(paths)
                            if ~isempty(paths{source})
                                sourceNumber = paths{source}(end);
                                len = numel(paths{source});
                                str = max(numel(paths{source})-obj.history,1);
                                historyToFeed = paths{source}(str:end);
                                transCost(source) = log(obj.NodeList{sourceNumber}.getTransitionProbToNode(targetNode,historyToFeed,{sourceNumber len}));
                            else
                                %transCost(source) = costPerPath(source); % question this
                                transCost(source) = log(obj.NodeList{source}.getTransitionProbToNode(targetNode,[],{target len}));
                            end
                        end
                        obj.TStore(target,:) = transCost(:);   
                    end
                    obj.needPRE = 0;
                end
                for target = 1:numel(obj.NodeList)
                    for source = 1:numel(paths)
                        [transitionCost(target) selPath(target)] = max(obj.TStore(target,:)+costPerPath);
                    end
                end
            end
            
            %transitionCost = bsxfun(@plus,transitionMatrix,costPerPath);
            %[transitionCost selPath] = max(transitionCost,[],2);
            
            
            for e = 1:numel(selPath)
                if isempty(paths{e})
                    newPaths{e} = [selPath(e) e];                    
                else
                    newPaths{e} = [paths{selPath(e)} e]; 
                end                
            end
            paths = newPaths;
            costPerPath = transitionCost;
        end
        
        function [num] = translateNodeToNumber(obj,node)
            num = [];
            for n = 1:numel(obj.NodeList)
                if obj.NodeList{n} == node
                    num = n;
                    break
                end
            end
        end
    end
    
end