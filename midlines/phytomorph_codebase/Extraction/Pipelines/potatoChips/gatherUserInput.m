function [] = gatherUserInput(I,rM)
    R = regionprops(logical(rM),'BoundingBox');
    for e = 1:numel(R)
        tmp = imcrop(I,R(e).BoundingBox);
        imshow(tmp,[]);
        button = questdlg('Stem End Defect?','Please Answer',default);
        
    end
end