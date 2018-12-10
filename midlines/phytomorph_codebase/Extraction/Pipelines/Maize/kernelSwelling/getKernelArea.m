function [boundary currentArea MASK] = getKernelArea(I,areaThresh,threshValues)   
    dispi = 0;
    try
        I = imfilter(I,fspecial('disk',7),'replicate');
        R = rgb2hsv_fast(I,'single','H');
        E = strel('disk',7);
        %MASK = R < 65/256 | R > 130/256;
        MASK = R < threshValues(1) | R > threshValues(2);
        MASK = bwareaopen(MASK,areaThresh);
        MASK = imclearborder(MASK);
        MASK = imclose(MASK,E);
        MASK = imfill(MASK,'holes');
        B = bwboundaries(MASK);
        boundary = B{1};
        currentArea = sum(MASK(:));
        if dispi
            imshow(MASK,[]);
            drawnow
        end
    catch ME
        boundary = [];
        currentArea = 0;
    end
end