function [HI CM BB] = sampleRing(rM,I)
    % sample the ring mask and return a histogram or pdf of colors
    % sample the symmetry of the color space.
    % non-symetrical color mean stem end defect
    R = regionprops(logical(rM),'PixelIdxList','Centroid','BoundingBox');
    
    for e = 1:numel(R)
        for c = 1:size(I,3)
            tmp = I(:,:,c);
            h = tmp(R(e).PixelIdxList);
            HI(:,c,e) = ksdensity(double(h),linspace(0,255,255));            
        end
        CM(:,e) = R(e).Centroid;
        BB(e,:) = R(e).BoundingBox;
    end
end