function [kidx] = gmmImage(I,gmm)
    I = permute(I,[3 1 2]);
    sz = size(I);
    I = reshape(I,[sz(1) prod(sz(2:3))])';
    kidx = gmm.cluster(I);
    kidx = reshape(kidx,sz(2:3));
end