function [Mout] = bwlarge(Min,n)
    if nargin == 1
        n = 1;
    end
    Mout = zeros(size(Min));
    R = regionprops(logical(Min),'PixelIdxList','Area');
    for e = 1:min(numel(R),n)
        R = regionprops(logical(Min),'PixelIdxList','Area');
        [J,midx] = max([R.Area]);
        
        Mout(R(midx).PixelIdxList) = 1;
        Min(R(midx).PixelIdxList) = 0;
    end
    Mout = logical(Mout);
end