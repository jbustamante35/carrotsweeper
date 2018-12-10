function [R] = getSingleMask(fileName)
    try
        I = imread(fileName);
        R = regionprops(logical(I),'Image');
    catch
        R = [];
    end
end