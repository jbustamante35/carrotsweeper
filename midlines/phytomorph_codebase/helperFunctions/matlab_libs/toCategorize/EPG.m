function [X] = EPG(P,X)
    %%% 
    % P := process graph store
    % X := data graph store
    %%%
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%
    % init eval of process graph
    %%% create node for query
    qn = node();
    qn.setProp('headProcess','1');    
    current_node = queryNodeProp(P,qn);
    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % while there is a next node in the process
    while (numel(current_node.l) == 2)        
         = hn.getData(X.getData);
        X.insertNode(nn);
        %%% assign the next node as current for iteration
        current_node = current_node.l{2}(2);
    end
        
    
    
end