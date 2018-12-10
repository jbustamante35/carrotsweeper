function [I] = rgbAndFilterExtract(I,filter)
    I = imfilter(I,filter,'replicate');

    I = permute(I,[3 1 2]);
    sz = size(I);
    I = reshape(I,[sz(1) prod(sz(2:3))])';
%{
    hsv = squeeze(rgb2hsv(permute(I/255,[1 3 2])));
    
    I = [I/255 hsv(:,2:3)];
    %}
end