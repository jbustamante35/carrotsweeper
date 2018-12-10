function [cropBOXES_Region,innerBlueSquare_Region] = getTypeCrops(I,map)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pad the edges with the red border
    for r = 1:4
        I(1,:,1) = 1;
        I(1,:,2) = 0;
        I(1,:,3) = 0;
        I = imrotate(I,90);
    end

    
    
    Lab = rgb2lab(I);
    
    %{
    % old method with map for cluster analysis
    % on the date of Dec.3, 2018 - recoding with Lab color space
    % and hard coding threshold for red-green > 25
    L_nodither = rgb2ind(I,map,'nodither');
    L_dither = rgb2ind(I,map,'dither');
    %RGB = label2rgb(L);
    redTape = L_dither == 4 | L_dither == 6 | L_dither == 2;
    redTape = bwlarge(redTape);
    redTape = imclose(redTape,strel('square',71));
    %}
    
    aChannel_threshold = 25;
    redTape = Lab(:,:,2) > aChannel_threshold;
    redTape = bwlarge(redTape);
    redTape = imclose(redTape,strel('square',71));
    
    
    cropBOXES = imfill(redTape,'holes') - redTape;
    cropBOXES = bwareaopen(cropBOXES,250000);
    cropBOXES = bwlarge(cropBOXES,12);
    cropBOXES_Region = regionprops(logical(cropBOXES),'BoundingBox','Area');

    
    
    cbArea = [];
    for i = 1:numel(cropBOXES_Region)
        cbArea(i) = cropBOXES_Region(i).BoundingBox(3)*cropBOXES_Region(i).BoundingBox(4);
        cbmaxDIM(i) = max(cropBOXES_Region(i).BoundingBox(3:4));
    end
    rm = cbArea > 700000 | cbmaxDIM > 800;
    cropBOXES_Region(rm) = [];



    %{ 
    % old method for blue channel - clustering on map
    % now threshold for Lab b channel < -30
    %}
    %blueSquare = L_dither == 7 | L_dither == 8;
    blueSquare = Lab(:,:,3) < -25;
    blueSquare = bwlarge(blueSquare);
    blueSquare = imclose(blueSquare,strel('square',71));
    blueSquare = imclearborder(blueSquare);


    innerBlueSquare = imfill(blueSquare,'holes') - blueSquare;
    innerBlueSquare = bwlarge(innerBlueSquare,1);
    innerBlueSquare = bwareaopen(innerBlueSquare,500000);
    innerBlueSquare_Region =  regionprops(logical(innerBlueSquare),'BoundingBox');
end