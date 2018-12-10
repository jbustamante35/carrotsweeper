function [OUT] = newSeedSizeSimple(I,bwareaValue,disp)

    % change to gray
    G = double(rgb2gray(I))/255;
    %thresh = graythresh(G);
    thresh = .5;
    B = G < thresh;
    
    B = bwareaopen(B,bwareaValue);
    
    R = regionprops(B,'Area','MajorAxis','MinorAxis','Image','Centroid','Orientation','PixelIdxList','Eccentricity','Perimeter','Extent','Image');
    
    dR = 1;
    
    A = [R.Area];
    [fidx1] = count(A);
    MA = [R.MajorAxisLength];
    [fidx2] = count(MA);
    mA = [R.MinorAxisLength];
    [fidx3] = count(mA);
    P = [R.Perimeter];
    [fidx4] = count(P);
    E = [R.Extent];
    [fidx5] = count(E);
    
    
    fidx = cat(1,fidx1,fidx2,fidx3,fidx4);
    fidx = fidx1;
    fidx = all(fidx,1);
    fidx = find(fidx);
    % seed single object image
    Z = zeros(size(B));
    for s = 1:numel(fidx)
        Z(R(fidx(s)).PixelIdxList) = 1;
    end
    out = flattenMaskOverlay(I,logical(B),[],'r');
    out = flattenMaskOverlay(out,logical(Z),.5,'g');
    
    if disp
        imshow(out,[]);
        waitforbuttonpress
    end
    
    % trim the regionprops
    R = R(fidx);
    
    
    F1 = [R.Area]*dR^-2;
    F2 = [R.MajorAxisLength]*dR^-1;
    F3 = [R.MinorAxisLength]*dR^-1;
    F4 = [R.Perimeter]*dR^-1;
    
    
    OUT = [F1' F2' F3' F4'];
end