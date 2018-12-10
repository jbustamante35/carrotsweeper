function [d] = curveFix(d)
    x = 1:size(d,1); 
    for e = 1:size(d,2)
        tmp = d(:,e);
        m = tmp == 0 | isinf(tmp);
        r = regionprops(m,'Area','PixelIdxList');
        fidx = find([r.Area]==1);
        for f = 1:numel(fidx)
            fsite = r(fidx(f)).PixelIdxList;
            tmpX = x;
            tmpX(fsite) = [];
            tmpY = tmp;
            tmpY(fsite) = [];
            tmp(fsite) = interp1(tmpX,tmpY,x(fsite));
        end
        d(:,e) = tmp;
    end
end