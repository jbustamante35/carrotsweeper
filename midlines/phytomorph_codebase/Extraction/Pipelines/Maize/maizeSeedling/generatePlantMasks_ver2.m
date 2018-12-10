function [mask] = generatePlantMasks_ver2(I,cluter_Level0,cluster_Level1,cluster_Level2)
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

    % this was working but i need to use the attached to base labels 
    % make the mask
    %mask = MT == 7 | MT == 5;
    % new make mask
    [mask] = getAttachedClusters(MT,[5 7],[4],[50 50]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % merge and cluster
    % find the mask
    qidx = find(mask);
    % cluster the mask
    [M_idxA] = cluster(cluster_Level2,[Q1(qidx,:)]);
    % convert mask to double
    mask = double(mask);
    % assign the clusters to the mask
    mask(qidx) = M_idxA;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find major objects the the obects next to them
    % find the base cluster type
    % changed to make cali work 1 or 2
    base = mask == 1 | mask == 2;
    % remove small objects
    base = bwareaopen(base,50);
    % find the
    bidx = find(base);
    % region props
    ridx = regionprops(imdilate(mask==2,strel('disk',1)),'PixelIdxList');
    % recolor near-by objects
    for e = 1:numel(ridx)
        if ~isempty(intersect(ridx(e).PixelIdxList,bidx))
            base(ridx(e).PixelIdxList) = 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % asign the output
    mask = base;
    mask = bwareaopen(mask,50);
    %mask = bwlarge(mask,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    
    
    
end