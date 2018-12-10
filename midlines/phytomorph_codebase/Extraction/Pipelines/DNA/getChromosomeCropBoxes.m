function [I,BOXES,rM,R] = getChromosomeCropBoxes(fileName)
    fprintf(['starting crop boxes.\n']);
    % read the image
    I = imread(fileName);
    if isdeployed()
        H = imfinfo(fileName);
        bitDepth = H.BitDepth/3;
    end
    % convert from 16 bit iamges
    I = double(I)*(2^bitDepth-1)^-1;
    % get the red channel
    R = I(:,:,1);
    % simple threshold on the red
    rM = R > graythresh(R);
    % remove small objects
    rM = bwareaopen(rM,100);
    % get the bounding boxes
    R = regionprops(rM,'BoundingBox','Centroid');
    % store boxes
    for e = 1:numel(R)
        BOXES{e} = R(e).BoundingBox;
    end
    % verbose
    fprintf(['ending crop boxes.\n']);
end