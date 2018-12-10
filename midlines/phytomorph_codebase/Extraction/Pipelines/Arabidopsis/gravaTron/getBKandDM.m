function [bk dM] = getBKandDM(I,f)
    for e = 1:size(I,4)
        G(:,:,e) = rgb2gray(double(I(:,:,:,e))/255);
        G(:,:,e) = imfilter(G(:,:,e),fspecial('disk',5),'replicate');
    end
    bk = mean(G(:,:,1:f),3);
    dM = bsxfun(@minus,G,bk);
end