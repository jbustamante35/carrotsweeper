function [newMask] = applyClusterStage(I,mask,gmm)

    fidx = find(mask);
    K = [];

    for k = 1:3
        tmp = I(:,:,k);
        K = [K tmp(fidx)];
    end

    kidx = cluster(gmm,K);


    newMask = double(mask);
    newMask(fidx) = kidx;


end