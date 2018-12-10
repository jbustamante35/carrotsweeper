function [preVARE_X,preVARE_Y,preVARE_LX,preVARE_LY,preX_model,norX_model,preY_model,X2Y] = pre_post_model(D,N,postDIM,preDIM)
    % split into pre and post treatement
    preData = D(1:N(1),:);
    postData = D(N(2):end,:);
    
    
    CORR = myHoldOut_pls(preData',postData',15,50,'HoldOut',.20);
    [v,lidx] = max(mean(CORR,2));
    %lidx = 12;
    [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(preData',postData',lidx);
    funcPre = @(X)([ones(size(X',1),1) X']*beta)';
    yPre = funcPre(preData);
    
    
    
    
    
    normalizedY = postData - yPre;
    
    [postS postC postU postE postL postERR postLAM] = PCA_FIT_FULL_T(normalizedY,size(normalizedY,1));
    preVARE_LY = postL(postDIM);
    
    
    [postS postC postU postE postL postERR postLAM] = PCA_FIT_FULL_T(normalizedY,postDIM);
    
    
    [preS preX preU preE preL preERR preLAM] = PCA_FIT_FULL_T(preData,size(preData,1));
    
    
    
    
    
    [wU,wS,wV] = svd(Weights);
    
    
    % not in pre space
    nF_WHOLE = bsxfun(@minus,preData,mean(preData,2));
    
    nF_PART = wU(:,1:lidx)*wU(:,1:lidx)'*nF_WHOLE;
    %tmpC = PCA_REPROJ_T(nF_WHOLE,wE,wU);
    %nF_PART = PCA_BKPROJ_T(tmpC,wE,wU);
    
    
    nF_LEFTOVER = nF_WHOLE - nF_PART;
    
    
    totalVarX = std(preData,1,2).^2;
    totalVarX = sum(totalVarX);
    
    
    totalVarY = std(postData,1,2).^2;
    totalVarY = sum(totalVarY);
    
    
    totalVar_predictedX = std(nF_PART,1,2).^2;
    totalVar_predictedX = sum(totalVar_predictedX);
    
    totalVar_predictedY = std(yPre,1,2).^2;
    totalVar_predictedY = sum(totalVar_predictedY);
    
    totalVar_unpredictedX = std(nF_LEFTOVER,1,2).^2;
    totalVar_unpredictedX = sum(totalVar_unpredictedX);
    
    totalVar_unpredictedY = std(normalizedY,1,2).^2;
    totalVar_unpredictedY = sum(totalVar_unpredictedY);
    
    
    preVARE_X = [totalVar_predictedX totalVar_unpredictedX] .* totalVarX.^-1;
    preVARE_Y = [totalVar_predictedY totalVar_unpredictedY] .* totalVarY.^-1;
    
    
    
    
    nF_PART = bsxfun(@plus,nF_PART,mean(preData,2));
    nF_LEFTOVER = bsxfun(@plus,nF_LEFTOVER,mean(preData,2));
    
    
    
    [nS nC nU nE nL nERR nLAM] = PCA_FIT_FULL_T(nF_LEFTOVER,size(nF_LEFTOVER,1));
    preVARE_LX(2) = nL(preDIM(2));
    [nS nC nU nE nL nERR nLAM] = PCA_FIT_FULL_T(nF_LEFTOVER,preDIM(1));
    
    
    [aS aC aU aE aL aERR aLAM] = PCA_FIT_FULL_T(nF_PART,size(nF_PART,1));
    preVARE_LX(1) = aL(preDIM(1));
    [aS aC aU aE aL aERR aLAM] = PCA_FIT_FULL_T(nF_PART,preDIM(2));
    
    
    
    preX_model.E = aE;
    preX_model.U = aU;
    preX_model.V = diag(aLAM);
    preX_model.func = @(X)(aE'*bsxfun(@minus,X,aU));
    preX_model.sim = @(X)(bsxfun(@plus,aE*X,aU));
    
    norX_model.E = nE;
    norX_model.U = nU;
    norX_model.V = diag(nLAM);
    norX_model.func = @(X)(nE'*bsxfun(@minus,X,nU));
    norX_model.sim = @(X)(bsxfun(@plus,nE*X,nU));
    
    preY_model.E = postE;
    preY_model.U = postU;
    preY_model.V = diag(postLAM);
    preY_model.func = @(preX,X)(postE'*bsxfun(@minus,X-funcPre(preX),postU));
    preY_model.sim = @(X)(bsxfun(@plus,postE*X,postU));
    
    X2Y = funcPre;
    close all
    
    
end