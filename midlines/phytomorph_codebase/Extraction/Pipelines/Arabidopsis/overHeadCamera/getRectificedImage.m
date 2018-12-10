function [nI,tform,resEstimate,tf] = getRectificedImage(I,GMModel,V,tform,resEstimate,disp)
    
    
    [tf,tform,resEstimate] = isCheckerBoard(I,tform,resEstimate,GMModel,V{5},500000,disp);
    


    W = round(.5*(size(I,1).*resEstimate^-1));
    NIP = round(2*W*resEstimate);
    [w1,w2] = ndgrid(linspace(-W,W,NIP),linspace(-W,W,NIP));
    XI = tform.transformPointsForward([w2(:),w1(:)]);
    nI = [];
    for k = 1:size(I,3)
        nI(:,k) = ba_interp2(double(I(:,:,k)),XI(:,1),XI(:,2));
    end
    nI = reshape(nI,[NIP NIP 3]);
    nI = nI / 255;
    for k = 1:3
        tmp = nI(:,:,k);
        m = min(tmp(:));
        tmp = tmp - m;
        M = max(tmp(:));
        tmp = tmp / M;
        nI(:,:,k) = tmp;
    end
end