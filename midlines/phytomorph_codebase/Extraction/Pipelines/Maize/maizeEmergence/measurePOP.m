function [mea fitSig] = measurePOP(rawSig,meaP)



    u = mean(rawSig,2);
    
    xlab = 1:size(u,1);
    [J xval] = min(abs(u - mean(u)));
    [para] = fminsearch(@(X)mySigmoid_ver0(xlab',X,u),[u(end) J xlab(xval)]); 
    [J,yp] = mySigmoid_ver0(xlab',para);
    
    g = gradient(yp);
    [mp tp] = max(g);
    WID = sum(g > mp*meaP);
    
    fitSig = [yp g];
    
    mea = [u(end) mp tp WID];
end