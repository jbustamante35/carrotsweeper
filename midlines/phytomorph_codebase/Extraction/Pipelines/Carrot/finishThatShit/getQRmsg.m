function [qrMSG] = getQRmsg(I,RED_SQUARE)

    for b = 1:numel(RED_SQUARE)
        MSK = (poly2mask(RED_SQUARE{b}(:,1),RED_SQUARE{b}(:,2),size(I,1),size(I,2)));
        R = regionprops(MSK);
        QR = imcrop(I,R(1).BoundingBox);
        qrMSG{b} = decode_qr(QR);
    end
end