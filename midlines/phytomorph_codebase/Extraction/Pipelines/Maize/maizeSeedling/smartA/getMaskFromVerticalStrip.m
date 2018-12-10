function [R] = getMaskFromVerticalStrip(vertStrip,imageSize,EXT,BKBOUND,I)
    % make an image mask
    MASK = repmat(vertStrip,[imageSize(1) 1]);
    % get the bounding boxes for each mask
    R = regionprops(MASK,'BoundingBox','PixelIdxList');
    for e = 1:numel(R)
        tmpMASK = zeros(size(MASK));
        tmpMASK(R(e).PixelIdxList) = 1;
        tmpMASK = imdilate(tmpMASK,strel('disk',EXT));
        tmpMASK(:,1:BKBOUND) = 0;
        tmpMASK(:,end-BKBOUND) = 0;
        tmpR = regionprops(logical(tmpMASK),'BoundingBox');
        R(e).BoundingBox = tmpR(1).BoundingBox;
    end
end