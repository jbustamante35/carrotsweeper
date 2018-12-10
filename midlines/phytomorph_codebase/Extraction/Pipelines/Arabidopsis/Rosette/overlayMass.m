function [] = overlayMass(I,BB,gmm,grp)
    I = imread(I);
    
    %% crop
    for e = 1:numel(BB)
        subI(:,:,:,e) = getSquarePots(I,BB{e});
        %ssubI(:,:,:,e) = imfilter(subI(:,:,:,e),fspecial('disk',3));
    end
    
    %% cluster pixels
    for e = 1:size(subI,4)
        cI = classifyImage(double(subI(:,:,:,e)),gmm);
        mask = cI == grp;
        %mask = logical(bwlarge(mask));
        mask = bwareaopen(mask,1000);
        mask = imclearborder(mask);
        out = flattenMaskOverlay(subI(:,:,:,e),mask);
        imshow(out,[]);
        waitforbuttonpress
    end
end