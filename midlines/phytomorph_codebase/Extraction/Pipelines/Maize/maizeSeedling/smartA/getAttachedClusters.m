function [base] = getAttachedClusters(maskIDX,baseIDX,attachedIDX,areaOpenValue)
    base = zeros(size(maskIDX));
    for e = 1:numel(baseIDX)
        base = base | maskIDX == baseIDX(e);
    end
    base = bwareaopen(base,areaOpenValue(1));
    attached = zeros(size(maskIDX));
    for e = 1:numel(attachedIDX)
        attached = attached | maskIDX == attachedIDX(e);
    end
    base = bwareaopen(base,areaOpenValue(2));
    % find the
    bidx = find(base);
    % region props
    ridx = regionprops(imdilate(attached,strel('disk',1)),'PixelIdxList');



    % recolor near-by objects
    for e = 1:numel(ridx)
        if ~isempty(intersect(ridx(e).PixelIdxList,bidx))
            base(ridx(e).PixelIdxList) = 1;
        end
    end

end