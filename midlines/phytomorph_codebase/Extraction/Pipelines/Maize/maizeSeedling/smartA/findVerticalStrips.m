function [eBLOCK] = findVerticalStrips(I)
    % hardcoded smooth value and top threshold
    smoothValue = 100;
    TOP_THRESH = 1150;
    % if image is char - read it
    if ischar(I)
        I = imread(I);
    end
    
    % crop off qr code
    I(1:TOP_THRESH,:,:) = [];
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the cone-tainers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: find the cone-tainer \n']);
    G = rgb2gray(single(I)/255);
    % smooth the image
    G = imfilter(G,fspecial('gaussian',[11 11]));
    % take grad
    [d1 d2] = gradient(G);
    % threshold grad along vertical direction
    d2 = abs(d2) > graythresh(abs(d2));
    %
    d2 = bwareaopen(d2,50);
    d2 = imclose(d2,strel('disk',10));
    sig = sum(abs(d2),1);
    sig = imfilter(sig,fspecial('average',[1 smoothValue]),'replicate');
    sig = bindVec(sig);
    threshSIG = graythresh(sig);
    % find the gaps
    BLOCK = sig > threshSIG;
    % remove the non-gaps that are less than 70
    BLOCK = bwareaopen(BLOCK,70);
    % close the pot holder chunks
    BLOCK = imclose(BLOCK,strel('disk',100));
    eBLOCK = BLOCK;
    %{
    eBLOCK = imdilate(BLOCK,strel('disk',EXT));
    % make an image mask
    MASK = ~repmat(eBLOCK,[size(I,1) 1]);
    % get the bounding boxes for each mask
    R = regionprops(~MASK,'BoundingBox','PixelIdxList');
    for e = 1:numel(R)
        tmpMASK = zeros(size(MASK));
        tmpMASK(R(e).PixelIdxList) = 1;
        tmpMASK = imdilate(tmpMASK,strel('disk',EXT));
        tmpMASK(:,1:BKBOUND) = 0;
        tmpMASK(:,end-BKBOUND) = 0;
        tmpR = regionprops(logical(tmpMASK),'BoundingBox');
        R(e).BoundingBox = tmpR(1).BoundingBox;
    end
    fprintf(['starting: find the cone-tainer :' num2str(numel(R)) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the cone-tainers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
end