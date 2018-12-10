function [startPoint] = findStartPoints(idxL,fileList)

    SNIP = 50;
    
    
    I = imread(fileList{1});
    imgSZ = size(I);
    Z = zeros(imgSZ);
    Z1 = Z;
    fidx1 = find(idxL(:,2) == 1);
    Z(idxL(:,1)) = 1;
    Z1(idxL(fidx1,1)) = 1;
    S1 = bwmorph(Z1,'skeleton',inf);
    
    
    
    
    U1 = mean(Z1,1);
    fidx = find(U1~=0);
    startColumn = fidx(end) - SNIP;
    Z(:,startColumn:end) = 0;
    S1(:,startColumn:end) = 0;
    S1 = bwareaopen(S1,100);
    
    
    R = regionprops(logical(S1(:,startColumn-1)),'PixelIdxList');
    R = regionprops(logical(S1),'PixelIdxList');
    
    for e = 1:numel(R)
        tmp = [];
        tmpMSK = zeros(size(I));
        tmpMSK(R(e).PixelIdxList) = 1;
        [tmp(:,1) tmp(:,2)] = find(tmpMSK);
        value = max(tmp(:,2));
        fidx = find(tmp(:,2)==value);
        startPoint(e,:) = mean(tmp(fidx,:),1);
    end
    
    
    
end