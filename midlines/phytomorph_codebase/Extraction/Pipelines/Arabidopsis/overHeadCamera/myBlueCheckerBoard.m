function [imagePoints,boardSize] = myBlueCheckerBoard(I)


    

    [imagePoints,boardSize] = detectCheckerboardPoints(rgb2gray(I));

    
    %{
    L_nodither = rgb2ind(I,map,'nodither');
    L_dither = rgb2ind(I,map,'dither');
    blueSquare = L_nodither == 7 | L_nodither == 8;
    %}
    
    blueSquare = bwlarge(blueSquare);
    
    innerBlueSquare = imfill(blueSquare,'holes') - blueSquare;
    innerBlueSquare = bwlarge(innerBlueSquare);
    R = regionprops(logical(innerBlueSquare),'BoundingBox');
    R = R.BoundingBox;
    POLY = [R(1:2);R(1:2) + [R(3) 0];R(1:2) + [R(3) R(4)];R(1:2) + [0 R(4)];R(1:2)];
    INNER = poly2mask(POLY(:,1),POLY(:,2),size(I,1),size(I,2));
    INNER = imerode(INNER,strel('square',41));

    
    Z = ba_interp2(double(INNER),imagePoints(:,1),imagePoints(:,2));
    imagePoints = imagePoints(Z==1,:);

    boardSize = [size(imagePoints,1)^.5 size(imagePoints,1)^.5];
    
    
end