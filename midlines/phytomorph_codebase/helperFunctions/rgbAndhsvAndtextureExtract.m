function [I] = rgbAndhsvAndtextureExtract(I,filter)
    % convert to double if possible
    if ~isa(I,'double')
        I = double(I);
    end

    G = rgb2gray(I/255);


    [d1 d2] = gradient(G);
    d = (d1.^2 + d2.^2).^.5;
    


    I = imfilter(I,filter,'replicate');

    I = permute(I,[3 1 2]);
    sz = size(I);
    I = reshape(I,[sz(1) prod(sz(2:3))])';

    hsv = squeeze(rgb2hsv(permute(I/255,[1 3 2])));
    
    I = [I/255 hsv(:,2:3) d(:)];
end