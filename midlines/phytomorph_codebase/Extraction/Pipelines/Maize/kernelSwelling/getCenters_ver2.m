function [centers] = getCenters_ver2(fileName,mainBOX)
    I = imread(fileName);
    I = imcrop(I,mainBOX);
    
    I = rgb2hsv(I);
    I = double(I(:,:,1));
    I = I < graythresh(I);
    I = imclearborder(I);
    I = bwareaopen(I,2000);
    S1 = std(I,1,1);
    bk1 = imerode(S1,ones(1,500));
    bk1 = imfilter(bk1,fspecial('disk',201));
    S1 = imfilter(S1,fspecial('disk',101));
    [lm1] = nonmaxsuppts(S1,601);
    R1 = regionprops(lm1,'Centroid');
    lm1 = zeros(size(lm1));
    for e = 1:numel(R1)
        lm1(round(R1(e).Centroid(1))) = 1;
    end
    S1 = bindVec(S1);
    level = graythresh(S1);
    lm1 = lm1 & S1 > level;
    
    S2 = std(I,1,2);
    bk2 = imerode(S2,ones(500,1));
    bk2 = imfilter(bk2,fspecial('disk',201));
    S2 = imfilter(S2,fspecial('disk',101));
    [lm2] = nonmaxsuppts(S2,601);
    R2 = regionprops(lm2,'Centroid');
    lm2 = zeros(size(lm2));
    for e = 1:numel(R2)
        lm2(round(R2(e).Centroid(2))) = 1;
    end
    S2 = bindVec(S2);
    level = graythresh(S2);
    lm2 = lm2 & S2 > level;
   
    M = double(lm2)*double(lm1);
    [c1 c2] = find(M);
    centers = [c1+mainBOX(2) c2+mainBOX(1)];
    
    
end