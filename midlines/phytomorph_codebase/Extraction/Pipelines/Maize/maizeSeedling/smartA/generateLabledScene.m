function [] = generateLabledScene(I,sceneStruct)


    if ischar(I)
        % read the image
        I = double(imread(I))/255;
    end
    
  
    [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(I,sceneStruct.bioDataBuffer);
    
    

    out = flattenMaskOverlay(bioStrip,bwlarge(sceneStruct.backgroundMask),.5,'r');
    out = flattenMaskOverlay(out,sceneStruct.coneTainerMask,.5,'b');

    
    
    outTOP = flattenMaskOverlay(dataStrip,logical(sceneStruct.qrMask),.5,'y');

    imshow([outTOP;out],[]);
    hold on
    for r = 1:numel(sceneStruct.conetainerBox)
        tmpBOX1 = sceneStruct.conetainerBox(r).BoundingBox;
        tmpBOX2 = sceneStruct.conetainerBox(r).MODBoundingBox;
        tmpBOX3 = sceneStruct.plantBox(r).BoundingBox;

        tmpBOX1(2) = tmpBOX1(2) + size(dataStrip,1);
        tmpBOX2(2) = tmpBOX2(2) + size(dataStrip,1);
        tmpBOX3(2) = tmpBOX3(2) + size(dataStrip,1);


        tmpBOX4 = sceneStruct.qrObject.cropBox;

        rectangle('Position',tmpBOX1,'EdgeColor','c');
        rectangle('Position',tmpBOX2,'EdgeColor','m');
        rectangle('Position',tmpBOX3,'EdgeColor','g');
        rectangle('Position',tmpBOX4,'EdgeColor','r');
    end
    drawnow
    hold off
end