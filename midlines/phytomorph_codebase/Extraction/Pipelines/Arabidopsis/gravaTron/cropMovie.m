function [cI] = cropMovie(I,BOX)
    cI = imcrop(I(:,:,:,end),BOX);
    cI = zeros([size(cI) size(I,4)]);
    for e = 1:size(I,4)
        cI(:,:,:,e) = imcrop(I(:,:,:,e),BOX);
    end
end