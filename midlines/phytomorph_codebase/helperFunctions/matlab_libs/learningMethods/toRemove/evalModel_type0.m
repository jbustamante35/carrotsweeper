function [learningVec] = evalModel_type0(model,data)
    
    % construct fibre space for data
    wholeCoeffs = PCA_REPROJ(data,model.basis,model.mean);
    wholeSim = PCA_BKPROJ(wholeCoeffs,model.basis,model.mean);
    wholeErr = sum((wholeSim - data).^2,2).^.5;
    wholeErr = whitenVec(wholeErr);    
    
    % construct learning vector
    wholeCoeffs = wholeCoeffs*model.lambda;
    learningVec = [wholeCoeffs,wholeErr];
    
    
end