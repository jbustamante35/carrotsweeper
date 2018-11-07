function [Z] = centroidFilter(B,thresh)
    CC = bwconncomp(B);
    R = regionprops(CC,'centroid','pixelidxList');
    Z = zeros(size(B));
    for e = 1:numel(R)
        cen = R(e).Centroid;
        if cen(2) > thresh & cen(1) > thresh
            Z(R(e).PixelIdxList) = 1;
        end
    end
end