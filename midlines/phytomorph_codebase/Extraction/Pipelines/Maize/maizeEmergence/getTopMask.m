function [nm] = getTopMask(i)
    u = rgb2hsv(i);
    nm = u(:,:,1) > .97 | u(:,:,1) < .065;
    nm = imfill(nm,'holes');
    nm = imopen(nm,strel('disk',11,0));
    nm = imclearborder(nm);
    nm = bwlarge(nm);
    nm = imfill(nm,'holes');
end