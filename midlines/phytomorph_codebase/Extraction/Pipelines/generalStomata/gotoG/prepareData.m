function [data] = prepareData(data,basisU,basisE)
    if ischar(data)
        [data] = myRemoteLoader(data,'T');
    end
    
    data.T = transForm_extractToFeatureTensor(data.T);
    
    
    
    
    
    
    
    
    
    
    
    
    
    data.C = [];
    for g = 1:3
        F = squeeze(data.T(:,:,g,:));
        szF = size(F);
        szF = [szF 1];
        F = reshape(F,[prod(szF(1:2)) szF(3)]);
        F = permute(F,[2 1]);
        data.C = [data.C PCA_REPROJ(F,basisE{g},basisU{g})];
    end
    data.C = data.C';
end