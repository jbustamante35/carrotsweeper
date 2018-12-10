function [skeletonIDX,absK] = findChromosomeEndPoints(contourPoints,skeletonPoints)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % boundary analysis to find the end points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the boundaries
    % filter the contour over scales
    fprintf(['Start:finding the skeleton end points.\n'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % measuring the curvature of the contour
    % finding peak curvature locations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshK = 100;
    % measure the curvature
    Kur = cwtK_closed_imfilter(contourPoints,{[41]});
    Kur = -(Kur.K);
    absK = sum(abs(Kur));
    kidx = (size(Kur,1)+1):(2*size(Kur,1));
    Kur = [Kur;Kur;Kur];
    dK = imdilate(Kur,strel('disk',51)) == Kur;
    dK = dK(kidx);
    cK = find(dK);
    Kur = Kur(kidx);
    ridx = abs(Kur(cK)) < threshK;
    cK(ridx) = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % snapping the curvature peaks to the nearest skeleton point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    skeletonIDX = [];
    for r = 1:numel(cK)
        skeletonIDX(r) = snapTo(skeletonPoints,contourPoints(cK(r),:));
    end
    fprintf(['End:finding the skeleton end points.\n'])
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OLD METHOD 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    % interesting method and might use this later
    % regions of higher curvature will move more under pressure of
    %smoothing
    S = 5:1:25;
    dBS = [];
    for N = 1:numel(S)
        dBS(:,:,N) = imfilter(contourPoints,ones(S(N),1)/S(N),'circular');
    end
    dF = diff(dBS,1,3);
    dF = mean(dF,3);
    %}