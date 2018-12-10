function [path] = traceSkeleton(skeletonImage,tipPoint)

    skeletonImage = bwmorph(skeletonImage,'skeleton',inf);
    ep = bwmorph(logical(skeletonImage),'endpoints',1);
    [e1 e2] = find(ep);
    PTS = [e1 e2];
    fidx = bsxfun(@eq,PTS,tipPoint);
    fidx = all(fidx,2);
    tipPoint = [e1(fidx) e2(fidx)];
    basePoint = [e1(~fidx) e2(~fidx)];


    [x,y] = find(skeletonImage);
    % stack the skeleton points for tracing
    DP = [x y]';
    % make adjaceny matrix
    T = Radjacency(DP,3);
    [idx(1)] = snapTo(DP',tipPoint);
    [idx(2)] = snapTo(DP',basePoint);
    [path , pathcost]  = dijkstra(T , idx(2) , idx(1));
    path = DP(:,path);
end