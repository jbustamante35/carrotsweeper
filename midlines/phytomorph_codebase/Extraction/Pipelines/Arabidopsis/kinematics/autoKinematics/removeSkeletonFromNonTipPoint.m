function [skeletonImage] = removeSkeletonFromNonTipPoint(skeletonImage,tipPoint,N)
    skeletonImage = bwmorph(skeletonImage,'thin',inf);
    skeletonImage = bwmorph(skeletonImage,'skeleton',inf);
    for k = 1:N
        ep = bwmorph(skeletonImage,'endpoints',1);
        [e1 e2] = find(ep);
        PTS = [e1 e2];
        fidx = bsxfun(@eq,PTS,tipPoint);
        fidx = all(fidx,2);
        tipPoint = [e1(fidx) e2(fidx)];
        e1(fidx) = [];
        e2(fidx) = [];
        skeletonImage(e1,e2) = 0;
    end
end