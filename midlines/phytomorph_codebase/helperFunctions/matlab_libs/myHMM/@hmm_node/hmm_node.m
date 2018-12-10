classdef hmm_node < handle

    properties
        StateName = '';
        
        NodeConnectionList = {};
        NodeTransitionProb = [];
        
        NodeReverseLookup = {};
        
        
        Distribution;
    end
    
    methods
        function [obj] = hmm_node(stateName)
            obj.StateName = stateName;
        end
        
        function [] = attachDistribution(obj,distribution,distributionNumber)
            obj.Distribution{distributionNumber} = distribution;
        end
        
        function [sample] = observe(obj,distributionNumber)
            sample = obj.Distribution{distributionNumber}.drawSample();
        end
        
        function [prob] = ProbObservation(obj,observation,distributionNumber,dn)
            prob = obj.Distribution{distributionNumber}.getProb(observation,dn);
        end
        
        function [] = attachNode(obj,node,prob)
            obj.NodeConnectionList{end+1} = node;
            %obj.NodeTransitionProb = [obj.NodeTransitionProb prob];
            obj.NodeTransitionProb{end+1} = prob;
        end
        
        function [node] = getNextNode(obj)
            idx = datasample(1:numel(obj.NodeConnectionList),1,'Weights',obj.NodeTransitionProb);
            node = obj.NodeConnectionList{idx};
        end
        
        function [prob] = forwardCalc(obj,NodeList,observation,observationLabel,dn)
            for e = 1:numel(NodeList)
                prob_emm(e) = NodeList{e}.getObservationProb(observation,observationLabel,dn);
                trans(e) = obj.toTransitionProb(NodeList{e});
            end
            prob = prob_emm.*trans;            
        end
        
        function [prob] = getObservationProb(obj,observation,observationLabel,dn)
            UQ = unique(observationLabel);
            prob = 1;
            for e = 1:numel(UQ)
                e;
                prob = prob * obj.ProbObservation(observation(observationLabel==UQ(e),:),UQ(e),dn(UQ(e)));
            end
        end
            
        function [prob] = toTransitionProb(obj,target)
            prob = 0;            
            for idx = 1:numel(obj.NodeConnectionList)
                if obj.NodeConnectionList{idx} == target
                    prob = obj.NodeTransitionProb(idx);
                    break                
                end
            end
        end
        
        function [] = upDateTransitionProb(obj,targetNode,transitionProb)
            for idx = 1:numel(obj.NodeConnectionList)
                if obj.NodeConnectionList{idx} == targetNode
                    obj.NodeTransitionProb{idx}.updateTransitionData(transitionProb);
                    break                
                end
            end
        end
        
        function [prob] = getTransitionProbToNode(obj,targetNode,stateObservation,augVars)
            idx = obj.getTargetNodeNumber(targetNode);
            if ~isempty(idx)
                prob = obj.NodeTransitionProb{idx}.computeTransitionProb(stateObservation,augVars);
            else
                prob = 0;
            end
        end
        
        function [idx] = getTargetNodeNumber(obj,targetNode)
            idx = [];
            for e = 1:numel(obj.NodeConnectionList)
                if obj.NodeConnectionList{e} == targetNode
                    idx = e;
                    break
                end
            end
        end
        
        
    end
    
    
end