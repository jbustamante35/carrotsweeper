function [qrData,boundingBox] = getCarrotQR(cellImage,kidx,qrClusterLabel)
    qrBOX = kidx == qrClusterLabel;
    qrBOX = bwlarge(qrBOX);
    qrBOX = imclose(qrBOX,strel('square',51));
    qrBOX = imfill(qrBOX,'holes') - qrBOX;
    qrBOX = bwlarge(qrBOX);
    R = regionprops(logical(qrBOX));
    boundingBox = R(1).BoundingBox;
    qrBOX = imcrop(cellImage,boundingBox);
    qrData = decode_qr(uint8(qrBOX));
end