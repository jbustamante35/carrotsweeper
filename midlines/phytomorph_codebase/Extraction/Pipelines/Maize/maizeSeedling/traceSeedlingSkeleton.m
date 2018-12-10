function [path,pathcost,midx,DP] = traceSeedlingSkeleton(skeleton,basePoint)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % trace the skeleton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: skeleton trace\n']);
    % find the skeleton for tracing
    [r c] = find(skeleton);
    % find the tips
    EP = imfilter(double(skeleton),ones(3,3));
    [re,ce] = find(EP == 2 & skeleton);
    % find skeleton
    [x,y] = find(skeleton);
    % stack the skeleton points for tracing
    DP = [x y]';
    fprintf(['starting: make adjacency matrix\n']);
    % make adjaceny matrix
    T = Radjacency(DP,(2^.5)+eps);
    fprintf(['ending: make adjacency matrix\n']);
    % find the longest path from the stem end point to the leaf tip
    pathcost = [];
    path = {};
    [idx(1)] = snapTo(DP',basePoint);
    for i = 1:numel(re)
        % find the end point in the skeleton
        [idx(2)] = snapTo(DP',[re(i) ce(i)]);
        fprintf(['starting: path trace\n']);
        % trace
        [path{i} , pathcost(i)]  = dijkstra(T , idx(1) , idx(2));
        fprintf(['ending: path trace\n']);
    end
    % set pathcost of inf to zero
    pathcost(isinf(pathcost)) = 0;
    % find the zeros - including inf path cost
    ridx = find(pathcost==0);
    % remove the 0 length paths
    pathcost(ridx) = [];
    % remove the 0 length paths
    path(ridx) = [];
    % find the max path cost
    [J,midx] = max(pathcost);
    fprintf(['ending: skeleton trace\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % trace the skeleton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end