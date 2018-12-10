function [rootCurve] = isolateRoot2(I,bulk)
    disp = 0;
    G = double(rgb2gray(I))/255;
    level = graythresh(G);
    B = G < level;
    R = regionprops(B,'Image','Area','BoundingBox','PixelIdxList');
    Area = [R.Area];
    [~,midx] = max(Area);
    B = zeros(size(B));
    B(R(midx).PixelIdxList) = 1;
    B = imfill(B,'holes');
    B = imdilate(B,strel('disk',bulk));
    P = B(1,:);
    fidx = find(B(1,:));
    dB = bwtraceboundary(B,[1 fidx(1)],'SW');
    dB(B(:,1) == 1,:) = [];
    rootCurve = dB;
end