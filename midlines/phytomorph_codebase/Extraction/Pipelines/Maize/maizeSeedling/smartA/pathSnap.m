function [value,sidx,path] = pathSnap(T,sourceIDX,targetsIDX)

    % find the longest path from the stem end point to the leaf tip
    pathcost = [];
    path = {};
    
    % trace from the source to the target(s)
    for i = 1:numel(targetsIDX)
        fprintf(['starting: path trace\n']);
        % trace
        [path{i} , pathcost(i)]  = dijkstra(T , sourceIDX , targetsIDX(i));
        fprintf(['ending: path trace\n']);
    end
    
   [value,sidx] =  min(pathcost(i));
   path = path{sidx};
    
    
    
end