function [midline,rootWidth] = getMidline(I)
    if ischar(I);I = imread(I);end
    skeletonImage = getSkeleton(I,0);
    tipPoint = findTipPoint(skeletonImage);
    skeletonImage = removeSkeletonFromNonTipPoint(logical(skeletonImage),tipPoint,100);
    path = traceSkeleton(skeletonImage,tipPoint);
    [midline,rootWidth] = extendMidline(I,path);
end