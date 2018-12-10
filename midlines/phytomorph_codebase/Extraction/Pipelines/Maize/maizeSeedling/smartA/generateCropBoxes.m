function [MASK] = generateCropBoxes(I)
   
    % find the vertical strips
    eBLOCK = findVerticalStrips(I);
    
    [R] = getMaskFromVerticalStrip(eBLOCK,size(I),50,20);
    [tmpD] = cropAndtrim(I,R,100);
    MASK = zeros(size(I,1),size(I,2));
    H = 3325 - 1000;
    for e = 1:numel(tmpD)
        
        tmpImage = imresize(tmpD{e}((end-H):end,:,:),[size(tmpD{e}((end-H):end,:,:),1) 1200]);
        
        
        
        Hi = double(tmpImage)/255;
        
        Hi = imfilter(Hi,fspecial('gaussian',[31 31],7),'replicate');
        E = edge(rgb2gray(Hi));
        E = imdilate(E,strel('disk',5,0));
        E(1:50,:) = 0;
        % find the 90-plumb
        [Hu,T,Ru] = hough(E','Theta',linspace(-5,5,100));
        P = houghpeaks(Hu,3);
        lines = houghlines(E',T,Ru,P,'FillGap',size(E,1)/2,'MinLength',500);
        xy = [lines(1).point1; lines(1).point2];
        v = round(mean(xy(:,1)));
        sig = zeros(size(Hi,1),1);
        sig(1:v) = 1;
        
        sig = repmat(sig,[1 size(tmpImage,2)]);
        sig = imresize(sig,[size(tmpD{e}((end-H):end,:,:),1) size(tmpD{e}((end-H):end,:,:),2)]);
        
        
        
        UL = round([R(e).BoundingBox(2) R(e).BoundingBox(1)]);
        
        MASK((end-H):end,UL(2):(UL(2)+size(sig,2)-1)) = sig;
       
    end
    MASK = sparse(logical(MASK));
end