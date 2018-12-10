function [MASK] = makeRootMask(I)
    % con
    tmp = rgb2hsv_fast(I,'single','S');

    
    h = fspecial('average',[5 5]);
    
    
    
    tmp = imfilter(tmp,h);
    t = graythresh(tmp);
    B = tmp < t;
    
    
    J = rgb2gray(double(I)/255);
    J = imfilter(J,h);
    delta = J - tmp;
    t = graythresh(delta);
    B = delta > t;
    
    
    OBJ = ~B;
    %OBJ = P(:,:,sel);
    OBJ = bwareaopen(OBJ,200);    
    
    
    R = regionprops(OBJ,'Area','PixelIdxList');
    
    
    R2 = regionprops(B,'Area','PixelIdxList');
    A = [R2.Area];
    [~,sidx] = max(A);
    
    % make B the background object - the largest object
    B = zeros(size(B));
    B(R2(sidx).PixelIdxList) = 1;
    
    
    fB = imfill(B,'holes');
    MASK = double(fB) - double(B);    
    MASK = logical(MASK);
    MASK = bwareaopen(MASK,10000);
    MASK = imopen(MASK,strel('disk',5));
    MASK = imdilate(MASK,strel('disk',101));
    
    rm = zeros(numel(R),1);
    for r = 1:numel(R)
        tmp = MASK(R(r).PixelIdxList);
        per = sum(tmp)/numel(R(r).PixelIdxList);
        if per < .70
            rm(r) = 1;
        end
    end
    R(find(rm)) = [];
    
    
    oMASK = zeros(size(MASK));
    for r = 1:numel(R)
        oMASK(R(r).PixelIdxList) = 1;
    end
        
    
    
    
    MASK = oMASK;
end