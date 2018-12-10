function [QRcode] = readQR(I,oPath,disp)

%black mask
QRMask = I(:,:,1) > 5 & I(:,:,1) < 105 & I(:,:,2) > 5 & I(:,:,2) < 105 & I(:,:,3) > 5 & I(:,:,3) < 105
R = regionprops(QRMask, 'boundingbox');
QRCode = imcrop(I, R(1).BoundingBox)
%read code

end
