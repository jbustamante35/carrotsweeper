function [qrLOCS_CORR_P,scaleOUT,sBOX,sSAMP] = fFind_QR(I,qrCenterNet,IMGSZ,sGRID,USE_SQUARE,ZPARA)

    sGRID = flip(sGRID,2);
    sBOX = point2Box(sGRID,USE_SQUARE);
    sSAMP = [];
    for e = 1:size(sBOX,1)
        sSAMP(:,:,:,e) = mimcrop(I,sBOX(e,:),IMGSZ);
    end
    qrLOCS = qrCenterNet.predict(sSAMP);
    qrLOCS = bsxfun(@times,qrLOCS,ZPARA(2,:));
    qrLOCS = bsxfun(@plus,qrLOCS,ZPARA(1,:));
    qrLOCS_CORR_P = [];
    for e = 1:size(qrLOCS,1)
        qrLOCS_CORR_P(e,:) = (qrLOCS(e,3))*((qrLOCS(e,1:2))) + (sBOX(e,1:2)) + .5*USE_SQUARE;
        scaleOUT(e) = qrLOCS(e,3);
    end
end
    
    