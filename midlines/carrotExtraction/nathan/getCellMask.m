function [cellMask,cellBoxes] = getCellMask(I,GMM,map,cellClusterNumber)
    smallThreshold = 500000;
    kidx = rgb2ind(I,map,'nodither');
    cellMask = 0;
    for e = 1:numel(cellClusterNumber)
        cellMask = cellMask | kidx == cellClusterNumber(e);
    end
    %cellMask = imclose(cellMask,strel('square',11));
    cellMask = bwareaopen(cellMask,smallThreshold);
    cellMask = imclearborder(cellMask);
    cellBoxes = regionprops(logical(cellMask),'BoundingBox');
end