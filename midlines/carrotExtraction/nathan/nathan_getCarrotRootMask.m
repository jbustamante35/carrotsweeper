function [carrotMask,carrotImage] = getCarrotRootMask(cellCrop,sub_kidx,cropLine)
    try
        % round one
        carrotImage = cellCrop(:,1:cropLine,:);
        carrotImage_kidx = sub_kidx(:,1:cropLine,:);
        HSV = rgb2hsv(carrotImage);
        carrotMask = HSV(:,:,2) > .2 & carrotImage_kidx ~= 2;
        carrotMask = bwlarge(carrotMask);
        iCM = sum(carrotMask,1);
        idx = find(iCM > 20);
        if ~isempty(idx)
            cropLine = idx(end);
            % round two to snap to edge
            carrotImage = cellCrop(:,1:cropLine,:);
            carrotImage_kidx = sub_kidx(:,1:cropLine,:);
            HSV = rgb2hsv(carrotImage);
            carrotMask = HSV(:,:,2) > .2 & carrotImage_kidx ~= 2;
            carrotMask = bwlarge(carrotMask);
        end
        midx = find(carrotMask(:,end));
        carrotMask(midx(1):midx(end),end) = 1;
        carrotMask = imfill(carrotMask,'holes');
    catch ME
        getReport(ME)
    end
end