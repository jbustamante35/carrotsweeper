function [m] = myStomataPeaks(m)

    md = bwdist(~m);
    
    
    
    
    for e = 1:100
        md = md - .01*del2(md);
    end
    
    
    m = imdilate(md,strel('disk',11,0)) == md;
    
    R = regionprops(m==1,'Area','Centroid','PixelIdxList');
    rm = [R.Area] ~= 1;
    R(rm) = [];
    
    
    m = zeros(size(m));
    m([R.PixelIdxList]) = 1;
    m = logical(m);
    midx = find(m);
    
    v = md(midx);
    
    ridx = v < 3;
    
    midx(ridx) = [];
    m = zeros(size(m));
    m(midx) = 1;
    m = logical(m);
end