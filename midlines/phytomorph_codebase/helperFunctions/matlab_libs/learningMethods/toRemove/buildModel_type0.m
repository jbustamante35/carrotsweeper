function [model] = buildModel_type0(data,para)
    %%%%%
    % construct fibre space for data
    [sim coeffs model.mean model.basis percentExp error model.lambda] = PCA_FIT_FULL(data,para.nComp); 
    model.org_lambda = model.lambda;
    model.lambda = diag(diag(model.lambda).^-.5);
end