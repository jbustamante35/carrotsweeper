function [BOX] = point2Box(pt,boxDim)
    for e = 1:size(pt,1)
        BOX(e,:) = [pt(e,:) - boxDim/2 boxDim];
    end
end