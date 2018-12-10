function [DP,T] = makeSkeletonTraceValues(skeleton)
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
end