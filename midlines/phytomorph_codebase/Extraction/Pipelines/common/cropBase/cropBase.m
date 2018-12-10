function [] = cropBase(image,U,E,G,Z,M)
    % image: image or image file
    % U: point in image to sample
    % E: reference frame to sample
    % G: grid to sample
    % Z: size of image after sample
    % M: method for sampling images
    
    if ischar(image)
        image = imread(image);
    end
    
    
end