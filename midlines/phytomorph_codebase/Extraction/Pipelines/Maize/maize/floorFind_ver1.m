function [BOX MASK tmp topFlag darkFlag] = floorFind_ver1(tmp)

    %% make background
    SNIP = 30;
    fSZ = 11;
    BK = [tmp(:,1:20) tmp(:,end-SNIP:end)];
    BK = mean(BK,2);
    BK = imfilter(BK,fspecial('average',[fSZ 1]),'replicate');
    BK = repmat(BK,[1 size(tmp,2)]);
    
    
    
    NHOOD = getnhood(strel('disk',5,0));
    sd = stdfilt(BK,NHOOD);
    sdValue = sum(sd(:));
    
    
    % mark top darkFlag
    topFlag = 0;
    if sdValue > (1.8*10^6)/1
        topFlag = 1;
    end
    
    
    
    
    
    %% darkFlag if dark and therefore need different thresholding
    darkFlag = 0;
    if mean(BK(:)) < 40
        tmp = log(double(tmp));
        darkFlag = 1;
    end
    darkFlag;
    tmp(isinf(tmp)) = 0;
    
    %% get the background again
    SNIP = 30;
    fSZ = 11;
    BK = [tmp(:,1:20) tmp(:,end-SNIP:end)];
    BK = mean(BK,2);
    BK = imfilter(BK,fspecial('average',[fSZ 1]),'replicate');
    BK = repmat(BK,[1 size(tmp,2)]);
    % line to remove the potential for 0 after log transform
    BK(isinf(BK)) = 0;
    
    
    
    
    
    %% make the kernel pop
    KERNEL = abs(double(tmp) - BK);
    KERNEL = bindVec(KERNEL);

    %% threshold if dark else hard code
    if darkFlag
        tv = linspace(0,1,200);
        for e = 1:numel(tv)
            MASK = KERNEL > tv(e);
            MASK = bwareaopen(MASK,100);
            tka(e) = sum(MASK(:));
        end
        %tka = log(tka);
        tka(isinf(tka)) = 0;
        ntka = tka/sum(tka);
        thresh = tv*ntka';
    else
        thresh = graythresh(KERNEL);
        thresh = .2;
    end
    
    
    %% make mask
    MASK = KERNEL > thresh;
    MASK = bwareaopen(MASK,200);
	MASK = imclose(MASK,strel('disk',10,0));
    MASK = imfill(MASK,'holes');
    
    
    %% if from the top - redo the mask
    topFlag;
    %
    if topFlag
        NHOOD = getnhood(strel('disk',21,0));
        %NHOOD = getnhood(strel('disk',5,0));
        toOP = edge(tmp,'canny');
        sd = stdfilt(toOP,NHOOD);
        M1 = sd < graythresh(sd);
        M2 = tmp > graythresh(tmp);
        MASK = M1.*M2;
        MASK = M1;
        MASK = bwareaopen(MASK,200);
        MASK = imfill(MASK,'holes');
        %MASK = imdilate(MASK,strel('disk',21,0));
        MASK = imdilate(MASK,strel('disk',5,0));
    end
    %
    
    hSTRIP = sum(MASK,2);
    hSTRIP = hSTRIP > 10;
    hSTRIP = bwareaopen(hSTRIP,50);
    fidx = find(hSTRIP);
    
    
    R = regionprops(MASK, 'BoundingBox','Area');
    [J,sidx] = max([R.Area]);
    %BOX = [0 0 size(tmp,2) fidx(end)];
    BOX = R(sidx).BoundingBox;
    
    
    %{
    out = flattenMaskOverlay(tmp, logical(MASK), .4);
    imshow(KERNEL,[]);
    stop = 1;
    %}
    
    
    
    
end