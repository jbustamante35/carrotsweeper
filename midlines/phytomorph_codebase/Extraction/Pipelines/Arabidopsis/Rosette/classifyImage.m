function [I] = classifyImage(I,gmm)
    sz = size(I);
    I = reshape(I,[prod(sz(1:2)) sz(3)]);
    I = cluster(gmm,I);
    I = reshape(I,sz(1:2));
end