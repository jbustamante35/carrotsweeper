function [nC] = getNewCenter(I,areaThresh)   
    dispi = 0;
    try
        
        I = imfilter(I,fspecial('disk',7),'replicate');
        R = rgb2hsv_fast(I,'single','H');
        
        
        
        % filter the background
        MASK = R < .2 | R > .5;
        %{
        R = rem(R + .4,1);
        level = graythresh(R);        
        MASK = R < mean(level);
        %}
        
        
        E = strel('disk',7);
        MASK = imclose(MASK,E);
        MASK = imfill(MASK,'holes');
        MASK = bwareaopen(MASK,areaThresh);
        R = regionprops(MASK,'Centroid');
        
        
        
        
        
        for e = 1:numel(R)
            delta(e) = norm(R(e).Centroid - [size(I,1) size(I,2)]/2);
        end
        [J,sidx] = min(delta);
        nC = R(sidx).Centroid - [size(I,1) size(I,2)]/2;
    catch ME
        getReport(ME)
        boundary = [];
        currentArea = 0;
    end
end