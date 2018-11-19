function [qrData,boundingBox] = getCarrotQR(cellImage,kidx,qrClusterLabel)
    qrBOX = kidx == qrClusterLabel;
    qrBOX = bwlarge(qrBOX);
    qrBOX = imfill(qrBOX,'holes') - qrBOX;
    qrBOX = bwlarge(qrBOX);
    R = regionprops(logical(qrBOX));
    boundingBox = R(1).BoundingBox;
    qrBOX = imcrop(cellImage,boundingBox);
    qrData = decode_qr(qrBOX);
end