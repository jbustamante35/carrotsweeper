function [I Z] = getChipMask_ver2(fileName,smallFill,chipFill)
% read image file
    I = imread(fileName);
    I = permute(I,[2 1 3]);
    G = rgb2gray(I);
    h = fspecial('gaussian',[11 11],5);
    Vo = imfilter(G,h,'replicate');
    level = 10;
    B = Vo > level;
    B = imclearborder(B);
    B = bwareaopen(B,smallFill);
    R = regionprops(B,'Area','PixelIdxList','Eccentricity');
    Z = zeros(size(B));
    
    
    E = [R.Eccentricity];
    fidx = find(E < .87);
    for r = 1:numel(fidx)
        Z(R(fidx(r)).PixelIdxList) = 1;
    end
    
    % close out small internal holes
    cZ = ~logical(Z);
    RzC = bwareaopen(cZ,chipFill);
    K = RzC ~= cZ;
    fidx = find(K);
    Z(fidx) = 1;
    
    %fix border
    Z = imerode(Z,strel('disk',5,0));
end
    