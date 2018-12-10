function [] = sorguhmLeafAnalysis(I,oPath)
    mkdir(oPath);
    [p,nm,ext] = fileparts(I);
    BUFFER_Y = 100;
    if ischar(I)
        I = imread(I);
    end
    
    G = rgb2gray(I);
    U2 = double(mean(G,2))/255;
    T = graythresh(U2);
    B = U2 > T;
    msk = repmat(B,[1 size(G,2)]);
    R = regionprops(logical(msk),'BoundingBox');
    IDX = reshape(1:numel(G),size(G));
    subIDX = imcrop(IDX,R(1).BoundingBox);
    subI = imcrop(I,R(1).BoundingBox);
    
    hsv = rgb2hsv_fast(subI);
    
    s1_1 = mean(hsv(:,:,1),1);
    b1 = s1_1 < graythresh(s1_1);
    s1_2= mean(hsv(:,:,3),1);
    b2 = s1_2 > graythresh(s1_2);
    msk2 = logical(b1.*b2);
    msk2 = repmat(msk2,[size(subI,1) 1]);
    msk2 = imerode(msk2,ones(1,40));
    leafMask = zeros(size(G));
    leafMask(subIDX(:)) = msk2(:);
    
    
    
    
    R1 = regionprops(logical(leafMask),'BoundingBox');
    R1(1).BoundingBox(2) = R1(1).BoundingBox(2) - BUFFER_Y;
    R1(1).BoundingBox(4) = R1(1).BoundingBox(4) + 2*BUFFER_Y;
    subI2 = imcrop(I,R1(1).BoundingBox);
    
    G2 = rgb2gray(double(subI2)/255);
    subM = G2 > graythresh(G2);
    HI = [];
    for k = 1:size(subI2,3)
        tmp = subI2(:,:,k);
        cl = double(tmp(find(subM)))/255;
        HI(k,:) = hist(cl,linspace(0,1,256));
    end
    widthProfile = sum(subM,1);
    
    G = rgb2gray(subI2);
    P = mean(G,2);
    mIDX = max(P);
    WID = 10;
    halfWidth = G(20:(mIDX-WID),:);
   
    
    image(I);
    hold on
    rectangle('Position',R1(1).BoundingBox,'EdgeColor','r');
    axis off
    X = round(R1(1).BoundingBox(1)):(round(R1(1).BoundingBox(1))+round(R1(1).BoundingBox(3)));
    Y = R1(1).BoundingBox(2)*ones(size(X)) - bindVec(widthProfile)*R1(1).BoundingBox(4);
    plot(X,Y);
    saveas(gca,[oPath nm '_output.tif']);
    csvwrite([oPath nm 'histogram.csv'],HI);
    close all
end

%{
    I = imread('/mnt/snapper/nate/forSorLeaf/SBFernandes_7-19-17/9K5B3062.CR2');
    sorguhmLeafAnalysis('/mnt/snapper/nate/forSorLeaf/SBFernandes_7-19-17/9K5B3062.CR2','./output/');
%}
