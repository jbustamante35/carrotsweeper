function [tipPoint] = findTipPoint(skeletonImage)
    Z = zeros(size(skeletonImage,1),size(skeletonImage,2));
    for r = 1:4
        Z(:,1) = 1;
        Z = imrotate(Z,90);
    end
    [z1 z2] = find(Z);
    tskel = logical(skeletonImage);
    tskel = bwmorph(tskel,'thin',inf);
    tskel = bwmorph(tskel,'skeleton',inf);
    
    
    ep = bwmorph(tskel,'endpoints',1);
    [e1 e2] = find(ep);
    for p = 1:2
        idx = snapTo([z1 z2],[e1(p) e2(p)]);
        dist(p) = (z1(idx) - e1(p)).^2 + (z2(idx) - e2(p)).^2;
    end
    
    [P,midx] = max(dist);
    tipPoint = [e1(midx) e2(midx)];
end