function [path] = traceMaizeMidline(skelPoints,startPoint,endPoint)
    % stack the skeleton points for tracing
    DP = skelPoints';
    fprintf(['starting: make adjacency matrix\n']);
    % make adjaceny matrix
    T = Radjacency(DP,3);
    fprintf(['ending: make adjacency matrix\n']);
    % find the longest path from the stem end point to the leaf tip
    pathcost = [];
    path = {};
    % snap the base point to the skeleton
    [idx(1)] = snapTo(DP',fliplr(endPoint));
    re = fliplr(startPoint);
    [idx(2)] = snapTo(DP',fliplr(startPoint));
    fprintf(['starting: path trace\n']);
    % trace
    [path , pathcost]  = dijkstra(T , idx(2) , idx(1));
    fprintf(['ending: path trace\n']);
    path = DP(:,path);
end