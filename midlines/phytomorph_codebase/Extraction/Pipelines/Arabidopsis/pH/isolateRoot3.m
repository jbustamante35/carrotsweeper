function [mask] = isolateRoot3(image,func,filtNumber,disp)
    
    eI = func(image);
    eI = bindVec(eI);    
    level = graythresh(eI);
    mask = eI < level;
    
    
    R = regionprops(mask,'PixelIdxList','Area');
    mA = [R.Area];
    [J,midx] = max(mA);
    mask = zeros(size(mask));
    mask(R(midx).PixelIdxList) = 1;
    mask = imdilate(mask,strel('disk',filtNumber-1));
    
    if disp
        imshow(mask);
        drawnow
    end
    
    
end