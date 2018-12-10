function [cI NL] = clusterImage(I,GMM)
    if ischar(I)
        I = double(imread(I))/255;
    end
    sz = size(I);
    I = reshape(I,[size(I,1)*size(I,2) size(I,3)]);
    [cI,NL] = cluster(GMM,I);
    cI = reshape(cI,sz(1:2));
end