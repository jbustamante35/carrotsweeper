function [I Z] = getChipMask_ver0(fileName,smallFill,chipFill)    
    I = imread(fileName);
    I = permute(I,[2 1 3]);
    %I = imresize(mag,.5);
    %{
    Vo = rgb2hsv_fast(I,'single','v');
    %}
    G = rgb2gray(I);
    h = fspecial('gaussian',[11 11],5);
    Vo = imfilter(G,h,'replicate');
    %BOX = 1.0e+04*[0.000651000000000   0.007851000000000   1.004398000000000   0.359998000000000];
    %[subI] = imcrop(V,BOX);
    %subI = V(100:end-100,round(1:size(V,2)*.4):end);
    %level = graythresh(subI);
    level = .35;
    level = 10;
    %H = hist(V(:),linspace(0,1,255));
    B = Vo > level;
    B = imclearborder(B);
    B = bwareaopen(B,smallFill);
    R = regionprops(B,'Area','PixelIdxList','Eccentricity');
    Z = zeros(size(B));
    
    
    E = [R.Eccentricity];
    fidx = find(E < .87);
    for r = 1:numel(fidx)
        Z(R(fidx(r)).PixelIdxList) = 1;
    end
    
    %{
    imshow(G,[]);
    hold on
    for e = 1:numel(R)
        GV(e) = mean(G(R(e).PixelIdxList));
    end
    
    
    
    
    fidx = find(GV < .7*255);
    fidx =
    for r = 1:numel(R)
        Z(R(fidx(r)).PixelIdxList) = 1;
    end
    %}
    
    
    
    %{
    E = [R.Eccentricity];
    fidx = find(E < .75);
    [fidx sel uF] = count(E);
    fidx = find(fidx);
    for r = 1:numel(fidx)
        Z(R(fidx(r)).PixelIdxList) = 1;
    end
    %}
    %{
    E = [R.Eccentricity];
    fidx = find(E < .5);
    for r = 1:numel(fidx)
        Z(R(fidx(r)).PixelIdxList) = 1;
    end
    %}
    %{
   
    %}
    
    
    % close out small internal holes
    cZ = ~logical(Z);
    RzC = bwareaopen(cZ,chipFill);
    K = RzC ~= cZ;
    fidx = find(K);
    Z(fidx) = 1;
    
    %fix border
    Z = imerode(Z,strel('disk',5,0));
end
    