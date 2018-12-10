function [level,G] = getThresholdLevel(I)
    % convert to hsv
    HSV = rgb2hsv(I);
    % get the gray scale image - the value component
    G = HSV(:,:,3);
    % smooth the image
    h = fspecial('gaussian',[31 31],5);    
    G = imfilter(G,h);
    % get the histogram
    [y xi] = hist(G(:),255);        
    % take the log of the histogram
    y = log(y);
    % dilate to find max peaks
    yd = imdilate(y,ones([1 31]));    
    fidx = find(yd == y);
    % get the height of the max peaks
    yvalue = y(fidx);
    % sort by max
    [~,sidx] = sort(yvalue,'descend');
    % get the location of the max
    fidx = fidx(sidx(1:2));
    % find the min between the peaks
    [~,loc3] = min(y(fidx(1):fidx(2)));
    % offset the min location by the first peak
    loc3 = loc3 + fidx(1);
    % get the level
    level = xi(loc3);
end