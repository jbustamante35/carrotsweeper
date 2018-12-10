function [MASK] = getObjectMask(I,areaThresh,toClear)
    I = imfilter(I,fspecial('disk',7));        
    cI = rgb2hsv(I);
    R = bindVec(double(cI(:,:,1)));
    level = graythresh(R);        
    MASK = R < mean(level);
    E = strel('disk',7);
    MASK = imclose(MASK,E);
    MASK = imfill(MASK,'holes');
    MASK = bwareaopen(MASK,areaThresh);
    if toClear
        MASK = imclearborder(MASK);
    end
end