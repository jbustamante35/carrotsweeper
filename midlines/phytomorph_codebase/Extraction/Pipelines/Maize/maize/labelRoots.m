function [gidx midx tipVec angle] = labelRoots(curve,pcaData,hmm,func,group,comp,disThresh)
    sz = size(curve.S);
    patchData = reshape(curve.S,[sz(1)*sz(2) sz(3)])';
    patchComp = PCA_REPROJ(patchData,pcaData.pE,pcaData.pU);
    TOT = [patchComp';bsxfun(@minus,curve.data,mean(curve.data,2))];
    observation_labels = [ones(size(patchComp,2),1);2*ones(2,1)];
    gidx = hmm.Viterbi(TOT,observation_labels,1);
    fidx = (gidx==group);
    [v midx] = func(patchComp(:,comp).*fidx');
    vec = bsxfun(@minus,curve.data,curve.data(:,midx));
    dis = sum(vec.*vec,1).^.5;
    sidx = dis < disThresh;
    [S C U tipVec L ERR LAM] = PCA_FIT_FULL(curve.data(:,sidx)',1);
    testVec = curve.data(:,midx) - U';
    if tipVec'*testVec < 0
        tipVec = -tipVec;
    end
    angle = atan2(tipVec(2),tipVec(1));
end