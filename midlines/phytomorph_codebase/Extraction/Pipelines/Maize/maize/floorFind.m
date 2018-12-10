function [BOX MASK tmp] = floorFind(tmp)

    %% make background
    SNIP = 30;
    fSZ = 11;
    BK = [tmp(:,1:20) tmp(:,end-SNIP:end)];
    BK = mean(BK,2);
    BK = imfilter(BK,fspecial('average',[fSZ 1]),'replicate');
    BK = repmat(BK,[1 size(tmp,2)]);
    
    %% flag if dark and therefore need different thresholding
    flag = 0;
    if mean(BK(:)) < 50
        tmp = log(double(tmp));
        flag = 1;
    end
    tmp(isinf(tmp)) = 0;
    
    %% get the background again
    SNIP = 30;
    fSZ = 11;
    BK = [tmp(:,1:20) tmp(:,end-SNIP:end)];
    BK = mean(BK,2);
    BK = imfilter(BK,fspecial('average',[fSZ 1]),'replicate');
    BK = repmat(BK,[1 size(tmp,2)]);
    BK(isinf(BK)) = 0;
    
    
    %% make the kernel pop
    KERNEL = abs(double(tmp) - BK);
    KERNEL = bindVec(KERNEL);

    %% threshold if dark else hard code
    if flag

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
        thresh = .3;
    end
    
    
    
    
    
    %% make mask
    MASK = KERNEL > thresh;
    MASK = bwareaopen(MASK,200);
	MASK = imclose(MASK,strel('disk',10,0));
    
    
    
    
    hSTRIP = sum(MASK,2);
    
    hSTRIP = hSTRIP > 10;
    hSTRIP = bwareaopen(hSTRIP,50);
    fidx = find(hSTRIP);
    BOX = [0 0 size(tmp,2) fidx(end)];
    
    MASK = imfill(MASK,'holes');
    
    %{
    out = flattenMaskOverlay(tmp, logical(MASK), .4);
    imshow(KERNEL,[]);
    stop = 1;
    %}
    
    
    
    %{
    E = edge(tmp,'canny');
    E = bwareaopen(E,10);
    [H,T,R] = hough(E,'Theta',-90);
    P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:)))) ;           
    lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
    max_len = 0;
    xy_long = [];
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];               

       %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       %plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       %plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    BOX = [0 0 size(tmp,2) xy_long(3)];
    %}
    
    
    
    
end