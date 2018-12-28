function [cellCrop,sub_kidx,qrData,QRboundingBox,cropLine] = getLeftRightBox(I,cellBoxes,map,gmm,blueClusterNumber,qrClusterLabel)


    for b = 1:numel(cellBoxes)
        cellCrop{b} = double(imcrop(I,cellBoxes(b).BoundingBox));
        [sub_kidx{b}] = gmmImage(cellCrop{b},gmm);
        [qrData{b},QRboundingBox{b}] = getCarrotQR(cellCrop{b},sub_kidx{b},qrClusterLabel);
        [cropLine(b)] = getCarrotSplitLine(sub_kidx{b},blueClusterNumber);
        
        
        if cropLine(b) > QRboundingBox{b}(1)
            cellCrop{b} = flip(cellCrop{b},2);
            sub_kidx{b} = flip(sub_kidx{b},2);
            [qrData{b},QRboundingBox{b}] = getCarrotQR(cellCrop{b},sub_kidx{b},qrClusterLabel);
            cropLine(b) = getCarrotSplitLine(sub_kidx{b},blueClusterNumber);
        end
        
        %{
        imshow(cellCrop{b}/255)
        drawnow
        %}
    end
end