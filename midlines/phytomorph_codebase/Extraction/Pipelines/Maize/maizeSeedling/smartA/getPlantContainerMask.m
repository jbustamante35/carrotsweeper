function [containerMask] = getPlantContainerMask(I,cluter_Level0,cluster_Level1)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create filtered image
    fI = imfilter(I,fspecial('average',[1 31]),'replicate');
    % create hsv image
    HSV = rgb2hsv(I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cast raw image - unfold tensor
    sz = size(I);
    Q1 = reshape(I,[prod(sz(1:2)) sz(3)]);
    % cast hsv image - unfold tensor
    sz = size(HSV);
    Q2 = reshape(HSV,[prod(sz(1:2)) sz(3)]);
    % cast filtered image - unfold tensor
    sz = size(fI);
    Q3 = reshape(fI,[prod(sz(1:2)) sz(3)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cluster at top level
    [M_idx,nlogl,P] = cluster(cluter_Level0,[Q1 Q2(:,[2 3])]);
    % resahpe the index clustering
    M_idx = reshape(M_idx,sz(1:2));
    % find the third cluster
    fidx = find(M_idx==3);
    % sub cluster the third group
    [M2_idx] = cluster(cluster_Level1,[Q3(fidx,:)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create master IDX image
    MT = M_idx;
    % assign the sub clusters
    MT(fidx) = M2_idx+3;
    containerMask = MT == 6;


















end