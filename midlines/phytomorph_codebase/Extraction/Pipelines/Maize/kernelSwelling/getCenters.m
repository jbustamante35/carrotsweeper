function [centers,uBBSZ,RxC,imageSZ,BOXSZ] = getCenters(fileName,mainBOX)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Start reading whole image\n']);
    I = imread(fileName);
    if nargin == 2
        I = imcrop(I,mainBOX);
    end
    fprintf(['End reading whole image\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the image size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imageSZ = size(I);
    % get the hue channel
    I = rgb2hsv_fast(I,'single','H');
    % filter the background
    I = I < .2 | I > .5;
    % remove small objects
    fprintf(['Start small object removal whole image\n']);
    I = bwareaopen(I,2000);
    fprintf(['End small object removal whole image\n']);
    fprintf(['Start clearing border\n']);
    I = imclearborder(I);
    fprintf(['End clearing border\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find obects via obect detection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Start object detection via count\n']);
    % get the regionprops
    R = regionprops(I,'Area','PixelIdxList','Perimeter','BoundingBox');
    % count objects on area
    [fidx1] = count([R.Area]);
    % count objects on perimeter
    [fidx2] = count([R.Perimeter]);
    % find those that match for both
    fidx = (fidx1==1) & (fidx2==1);
    fidx = find(fidx);
    R = R(fidx);
    for r = 1:numel(R)
        BOXSZ(r,:) = R(r).BoundingBox(3:4);
    end
    BOXSZ = max(BOXSZ,[],1);
    BOXSZ = round((BOXSZ + .5*BOXSZ)/2);
    BOXSZ = [BOXSZ(1) BOXSZ(1) BOXSZ(2) BOXSZ(2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MASK = zeros(size(I));
    for e = 1:numel(fidx)
        MASK(R(e).PixelIdxList) = 1;
        BB(e,:) = R(e).BoundingBox;
    end
    uBBSZ = mean(BB(:,3:4));
    fprintf(['End object detection via count\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find peaks along the horizontal 
    % std is large via graythresh and peaks via nonmax
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S1 = std(MASK,1,1);
    %bk1 = imerode(S1,ones(1,500));
    %bk1 = imfilter(bk1,fspecial('disk',201),'replicate');
    S1 = imfilter(S1,fspecial('disk',101),'replicate');
    [lm1] = nonmaxsuppts(S1,100);
    R1 = regionprops(lm1,'Centroid');
    lm1 = zeros(size(lm1));
    for e = 1:numel(R1)
        lm1(round(R1(e).Centroid(1))) = 1;
    end
    S1 = bindVec(S1);
    level = graythresh(S1);
    lm1 = lm1 & S1 > level;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find peaks along the vertical 
    % std is large via graythresh and peaks via nonmax
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S2 = std(MASK,1,2);
    %bk2 = imerode(S2,ones(500,1));
    %bk2 = imfilter(bk2,fspecial('disk',201),'replicate');
    S2 = imfilter(S2,fspecial('disk',101),'replicate');
    [lm2] = nonmaxsuppts(S2,101);
    R2 = regionprops(lm2,'Centroid');
    lm2 = zeros(size(lm2));
    for e = 1:numel(R2)
        lm2(round(R2(e).Centroid(2))) = 1;
    end
    S2 = bindVec(S2);
    level = graythresh(S2);
    lm2 = lm2 & S2 > level;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find centers in assumed square
    % count rows and columns
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = double(lm2)*double(lm1);
    [c1,c2] = find(M);
    centers = [c1 c2];
    numRows = sum(lm2);
    numCols = sum(lm1);
    RxC = [numRows numCols];
    %centers = [c1+mainBOX(2) c2+mainBOX(1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%{


%{
    %I = imread(fileName);
    %{
    %I = imcrop(I,mainBOX);
    %I = rgb2hsv(I);
    I = rgb2hsv_fast(I,'single','H');
    
    %{
    % simpler for problem
    fI = imfilter(I,fspecial('disk',31),'replicate');
    BK = imdilate(imresize(fI,.25),strel('disk',61,0));
    BK = imfilter(BK,fspecial('average',211),'replicate');
    BK = imresize(BK,size(I));
    I = I - BK;
    I = bindVec(I);
    % simpler for problem
    %}
    
    %I = double(I(:,:,1));
    
    %}
%}
    
%}