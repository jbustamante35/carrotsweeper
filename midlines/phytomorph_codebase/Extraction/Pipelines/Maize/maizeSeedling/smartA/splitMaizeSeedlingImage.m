function [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(I,bufferZone)
    [msg,qrCropBox] = getQRcode(I);
    cropLine = round(qrCropBox(2) + qrCropBox(4) + bufferZone);
    dataStrip = I(1:cropLine,:,:);
    bioStrip = I((cropLine+1):end,:,:);
end