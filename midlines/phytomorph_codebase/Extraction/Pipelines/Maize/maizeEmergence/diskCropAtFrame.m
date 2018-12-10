function [frameCrop] = diskCropAtFrame(imageStack,frameNumber,cropBox,BORDER)
    cropBox(1:2) = cropBox(1:2) - BORDER;
    cropBox(3:4) = cropBox(3:4) + 2*BORDER;
    cropBox = round(cropBox);
    for e = 1:numel(frameNumber)
        frameCrop(:,:,:,e) = imread(imageStack{frameNumber(e)},'PixelRegion',{[cropBox(2) cropBox(2)+cropBox(4)],[cropBox(1) cropBox(1)+cropBox(3)]});
    end
end