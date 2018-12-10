function [] = displayMaskNumbers(MASK,imageFile)
    tmpI = imread(imageFile);
    R = regionprops(logical(MASK),'Centroid');
    out = flattenMaskOverlay(tmpI, logical(MASK),.15,'g');
    imshow(out,[]);
    hold on
    for e = 1:numel(R)
        text(R(e).Centroid(1),R(e).Centroid(2),num2str(e));
    end
end