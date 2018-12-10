function [] = zoomGather(T,BOX)
    while T.hasdata
        [I,imageKey] = T.read();
        [xy(:,1),xy(:,2),V] = impixel(I);
        for e = 1:size(xy,1)
            b = point2Box(xy(e,:),BOX);
            subI = imcrop(I,b);
            imshow(subI,[]);
            drawnow
        end
    end
    
end