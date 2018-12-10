function [I,subI,boundingBoxes,connectedMASK,MASK,SKELETON] = cropImages_ver2(I,TOP_THRESH,cluter_Level0,cluster_Level1,cluster_Level2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read image if file name passed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(I)
        tic;
        fprintf(['start::reading file \n']);
        I = double(imread(I));
        [I angle] = rectifyImage(I/255);
        fprintf(['end::reading file::' num2str(toc) '\n']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read image if file name passed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SNAP = 50;
    I(:,(end-SNAP:end),:) = [];

    
    bottomImage = I(TOP_THRESH:end,:,:);
    [mask] = generatePlantMasks_ver2(bottomImage,cluter_Level0,cluster_Level1,cluster_Level2);
    
    R = regionprops(mask,'BoundingBox');
    
    for e = 1:numel(R)
        boundingBoxes{e} = R(e).BoundingBox;
        boundingBoxes{e}(2) = boundingBoxes{e}(2) + TOP_THRESH;
        subI{e} = imcrop(bottomImage,R(e).BoundingBox);
        MASK{e} = imcrop(mask,R(e).BoundingBox);
        connectedMASK{e} = connectPlant(MASK{e});
        SKELETON{e} = getPlantSkelton(MASK{e});
    end
    
    
    %{
    
     tic;
        fprintf(['starting:: get mask \n']);
        % get the mask
        MASK{e} = getMASK_ver0(I{e});
        % connect the mask
        connectedMASK{e} = connectPlant(MASK{e});
        % get the largest object in the mask
        connectedMASK{e} = bwlarge(connectedMASK{e});
        fprintf(['ending:: get mask::' num2str(toc) '\n']);
        tic
        fprintf(['starting:: get skeleton \n']);
        SKELETON{e} = getPlantSkelton(MASK{e});
        fprintf(['ending:: get skeleton::' num2str(toc) '\n']);
%}
end










