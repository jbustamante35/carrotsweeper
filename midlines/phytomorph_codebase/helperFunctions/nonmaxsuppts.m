function [localMAX] = nonmaxsuppts(I, rad)
    %%%%%%%%%%%%%%%%%%%%%
    % find the local corner
    %%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I           : = image for finding the local max    
    %           rad         : = the size of the area for local max
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           localMAX    : = binary map with local max values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dilate the image
    mx = imdilate(I,ones(rad));
    % local max map
    localMAX = (I == mx);
    %{
    %%%%%%%%%%
    % find all local max points
    fidx = find(localMAX);
    % sample the corner values at the local max values
    mI = cim(fidx);
    [nmI para] = normalizeMap(mI);
    level = graythresh(nmI);
    mI = normalizeMap(cim,para);
    clusterI = mI > level;
    % Find maxima, threshold, and apply bordermask
    cimmx = localMAX & clusterI & (cim>thresh) & MSK;
    %}
end