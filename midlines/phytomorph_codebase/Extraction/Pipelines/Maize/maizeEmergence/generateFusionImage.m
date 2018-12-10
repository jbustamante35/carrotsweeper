function [fusionImage] = generateFusionImage(rawImage,modelImage,clusterTree,mSelect)
    k = clusterTree.clusterImage(modelImage);
    % create mask
    mask = zeros(size(k));
    for s = 1:numel(mSelect)
        mask = mask | k == mSelect(s);
    end
    %mask = imclearborder(mask);
    mask = bwlarge(mask);
    % find fusion locations
    fidx = find(mask);
    fusionImage = modelImage;
    for e = 1:size(modelImage,3)
        tmpM = modelImage(:,:,e);
        tmpR = rawImage(:,:,e);
        tmpM(fidx) = tmpR(fidx);
        fusionImage(:,:,e) = tmpM;
    end
end