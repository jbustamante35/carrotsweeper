function [C] = cropCells(I,R)
    for e = 1:numel(R)
        C{e} = imcrop(I,R(e).BoundingBox);
    end
end