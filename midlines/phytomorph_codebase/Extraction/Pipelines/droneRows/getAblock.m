function [block] = getAblock(I,BOX)
    block= imcrop(I,BOX);
end