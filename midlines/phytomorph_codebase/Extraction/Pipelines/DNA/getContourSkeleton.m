function [contourPoints,skeletonPoints] = getContourSkeleton(oMASK,padValue,scaleFactor,meltIter,alpha)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % crop distance transform for contour and midline interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sz = size(oMASK);
    % pad with zeros
    PD = padarray(oMASK,padValue,0,'both');
    % generate distance transform
    DIST = bwdist(~PD) - bwdist(PD);
    % resize the distance transform
    DIST = imresize(DIST,scaleFactor);
    % melt the distance transform
    for t = 1:meltIter
        DIST = DIST + alpha*del2(DIST);
    end
    % get the contour from the melted transform
    contourPoints = contourc(double(DIST),[0 0]);
    % get the contour from the contour matrix
    contourPoints = contourPoints(:,2:contourPoints(2,1))';
    % make the new mask
    newSCM = poly2mask(contourPoints(:,1),contourPoints(:,2),(sz(1)+2*padValue(1))*scaleFactor,(sz(2)+2*padValue(2))*scaleFactor);
    % skeleton
    skeletonS = bwmorph(newSCM,'thin',inf);
    % find the skelton
    [skeletonPoints(:,2),skeletonPoints(:,1)] = find(skeletonS);
    skeletonPoints = bsxfun(@minus,(skeletonPoints * scaleFactor^-1),padValue);
    % scale and translate the contour
    contourPoints = bsxfun(@minus,(contourPoints*scaleFactor^-1),padValue);
    %{
    % display after contour and midline
    figure;
    imshow(tmpI,[]);
    hold on
    plot(dC(:,1),dC(:,2),'y')
    plot(skeletonPoints(:,1),skeletonPoints(:,2),'y.')
    hold off
    %waitforbuttonpress
    %}
end