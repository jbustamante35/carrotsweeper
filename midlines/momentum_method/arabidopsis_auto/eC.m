function [] = eC(imageName)
    mag = 2;
    I = double(imread(imageName))/255;
    I(1,:) = I(2,:);
    Io = I;
    I = cat(2,fliplr(I(:,1:SNIP)),I);
    I = imfilter(I,fspecial('disk',11),'replicate');
    %[g1 g2] = gradient(I);   
    %gV = (g1.^2 + g2.^2).^.5;
    I = abs(I-1);
    gV = I;
    gV = bindVec(gV);
    gV = imresize(gV,1/mag);
    for i = 1:400
        gV = gV + mydel2(gV);
        gV = bindVec(gV);
        imshow(gV,[])
        drawnow
        i
    end
    gV = imresize(gV,size(I));
    I = gV(:,SNIP+1:end);
    I = bindVec(I);
    gV = bindVec(gV);
    R = Io+Io.*gV;       
    gradThresh = graythresh(gV);    
    BLOB = gV > gradThresh;
    
     
    R = regionprops(BLOB,'Area','PixelIdxList');

    [J sidx] = sort([R.Area]);
    BLOB = zeros(size(I));                
    BLOB(R(sidx(end)).PixelIdxList) = 1;
    
    
    skel = bwmorph(BLOB,'skel',inf);
    skel = bwmorph(skel,'spur',inf);
    skel = imdilate(skel,strel('disk',21));
    B = bwboundaries(skel);
    imshow(cat(3,Io,Io,Io));
    hold on
    plot(B{1}(:,2),B{1}(:,1),'r')
    hold off
    drawnow
    saveas(gca,fn);
    
end