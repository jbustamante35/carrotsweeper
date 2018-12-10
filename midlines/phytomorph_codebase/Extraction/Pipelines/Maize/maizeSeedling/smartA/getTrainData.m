function [MASK,sI,sM,rI] = getTrainData(fileName,newY)
    % read images
    I = imread(fileName);
    [I angle] = rectifyImage(I);
    % get masks
    MASK = generateCropBoxes(I);
    % crop 
    [sI,sM] = verticalCrop_ver0(I,full(MASK),newY);
    % re-size
    rI = imresize(I,.25);
end