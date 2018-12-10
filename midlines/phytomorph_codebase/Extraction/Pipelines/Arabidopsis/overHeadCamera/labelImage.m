function [label] = labelImage(I,GMModel)
    sz = size(I);
    hsv = rgb2hsv(double(I)/255);
    I = reshape(I,[prod(sz(1:2)) sz(3)]);

    %hsv = rgb2hsv(double(I)/255);
    %I = double(I);
    %[label,~,prob] = GMModel.cluster([I(:,2)]);

    %label = reshape(label,sz(1:2));

    
    
    
    %greenPlant = label==4 & hsv(:,:,1) > .15 & hsv(:,:,1) < .3 & hsv(:,:,2) > .3;% & hsv(:,:,2) < .8 
    greenPlant = hsv(:,:,1) > .15 & hsv(:,:,1) < .388 & hsv(:,:,2) > .3;% & hsv(:,:,2) < .8
    greenPlant = imclose(greenPlant,strel('disk',3,0));
    greenPlant = bwareaopen(greenPlant,20);
    
    redTape = (hsv(:,:,1) < .05 | hsv(:,:,1) > .97) & hsv(:,:,3) > .4;
    %greenPlant = (hsv(:,:,1) > .25 & hsv(:,:,1) < .41) & hsv(:,:,2) > .189;
    whiteObjects = hsv(:,:,2) < .3 & hsv(:,:,3) > .8;
    blueBorder = (hsv(:,:,1) > .3 & hsv(:,:,1) < .7) & hsv(:,:,2) > .38 & hsv(:,:,3) < .77;
    blackBackground = ~redTape & ~greenPlant & ~whiteObjects & ~blueBorder;
    label = zeros(prod(sz(1:2)),1);
    label(find(blackBackground)) = 5;
    label(find(redTape)) = 1;
    label(find(whiteObjects)) = 3;
    label(find(blueBorder)) = 4;
    label(find(greenPlant)) = 2;
    %label = vec2ind([redTape(:) greenPlant(:) whiteObjects(:) blueBorder(:) blackBackground(:)]');
    
    label = reshape(label,sz(1:2));

    %probBorder = reshape(prob(:,4),sz(1:2));
end