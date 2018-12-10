function [sm,mm] = squareMini(s,t)
    
    % make new mask
    u = mean(s(:,:,:,1:t),4);
    u = rgb2hsv(u);
    nm = u(:,:,1) > .97 | u(:,:,1) < .06;
    %nm = imopen(nm,strel('disk',11,0));
    nm = imclearborder(nm);
    nm = bwlarge(nm);
    nm = imfill(nm,'holes');
    
    
    
    R = regionprops(nm,'BoundingBox','EquivDiameter');
    d = round(0.1272*R(1).EquivDiameter);
    
    
    m = imdilate(logical(nm),strel('disk',d,0));
    sm = zeros(200,200,3,size(s,4));
    R = regionprops(logical(m),'BoundingBox','EquivDiameter');
    
    for e = 1:size(s,4)
        tmp = imcrop(s(:,:,:,e),R(1).BoundingBox);
        sm(:,:,:,e) = imresize(tmp,[200 200])/255;
    end
    
   
    tmp = imcrop(m,R(1).BoundingBox);
    mm = imresize(tmp,[200 200]);
end