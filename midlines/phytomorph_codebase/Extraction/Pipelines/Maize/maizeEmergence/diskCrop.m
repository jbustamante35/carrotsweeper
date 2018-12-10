function [miniStack cropMask] = diskCrop(imageStack,MASK,cropBox,BORDER,rec)
    % init miniStack
    miniStack = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % expand each crop box by border amount
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(cropBox)
        cropBox{e}(1:2) = cropBox{e}(1:2) - BORDER;
        cropBox{e}(3:4) = cropBox{e}(3:4) + 2*BORDER;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % crop the mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    reSZ = [200 200];
    for e = 1:numel(cropBox)
        % crop the mask
        tmp = imcrop(MASK,cropBox{e});
        % transform the mask to overhead perspective
        tmp = imtransform(tmp,rec,'UData',[cropBox{e}(1) cropBox{e}(1)+cropBox{e}(3)],'VData',[cropBox{e}(2) cropBox{e}(2)+cropBox{e}(4)]);
        % clear the border
        tmp = imclearborder(tmp);
        % get the largest object
        tmp = bwlarge(tmp);
        % fill in holes
        tmp = imfill(tmp,'holes');
        % get the bounding box and equivdiameter
        R = regionprops(tmp,'BoundingBox','EquivDiameter');
        % scale for dilate
        d = round(0.1272*R(1).EquivDiameter);
        % dilate the mack
        tmp = imdilate(logical(tmp),strel('disk',d,0));
        % get the bounding box
        R = regionprops(logical(tmp),'BoundingBox','EquivDiameter');
        % sub crop
        tmp = imcrop(tmp,R(1).BoundingBox);
        
        tmp = imresize(tmp,reSZ);
        
        cropMask(:,:,e) = tmp;
        cropBox{e} = round(cropBox{e});
        subcropBox{e} = R(1).BoundingBox;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test crop for pre allocate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %testC = imread(imageStack{1},'PixelRegion',{[cropBox(2) cropBox(2)+cropBox(4)],[cropBox(1) cropBox(1)+cropBox(3)]});
    %testC = imtransform(testC,rec,'UData',[cropBox(1) cropBox(1)+cropBox(3)],'VData',[cropBox(2) cropBox(2)+cropBox(4)]);
    miniStack = zeros([reSZ 3 (numel(imageStack)-1) numel(cropBox)]);
    
    for e = 2:numel(imageStack)
        reportLoop('Cropping images',e-1,30);
        I = imread(imageStack{e});
        % loop over crop box and crop images
        for c = 1:numel(cropBox)
             tmp = imcrop(I,cropBox{c});
             tmp = imtransform(tmp,rec,'UData',[cropBox{c}(1) cropBox{c}(1)+cropBox{c}(3)],'VData',[cropBox{c}(2) cropBox{c}(2)+cropBox{c}(4)]);
             tmp = imcrop(tmp,subcropBox{c});
             tmp = imresize(tmp,reSZ);
             miniStack(:,:,:,e-1,c) = tmp;
        end
        
        
        
        %tmp = imread(imageStack{e},'PixelRegion',{[cropBox(2) cropBox(2)+cropBox(4)],[cropBox(1) cropBox(1)+cropBox(3)]});
        
        
        %miniStack(:,:,:,e-1) = imtransform(tmp,rec,'UData',[cropBox(1) cropBox(1)+cropBox(3)],'VData',[cropBox(2) cropBox(2)+cropBox(4)]);
    end
    
    %cropMask = imtransform(cropMask,rec,'UData',[cropBox(1) cropBox(1)+cropBox(3)],'VData',[cropBox(2) cropBox(2)+cropBox(4)]);
end