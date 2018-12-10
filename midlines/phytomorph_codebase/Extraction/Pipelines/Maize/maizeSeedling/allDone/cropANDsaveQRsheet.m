function [] = cropANDsaveQRsheet(I,pointList,BUFFER,oPath,oName)
    W = pointList(3,1) - pointList(1,1);
    H = pointList(7,2) - pointList(1,2);
    qrSheet = imcrop(I,[pointList(1,:)-BUFFER [W H] + 2*BUFFER]);
    
    imwrite(qrSheet,[oPath oName '.tif']);

    imshow(qrSheet,[]);
    drawnow


end