function [cData,U,E] = ncPCA(D,N)
    cData = [];
    for e = 1:numel(D)
        [S C U{e} E{e} L ERR LAM] = PCA_FIT_FULL(D{e},min(N(e),size(D{e},2)));
        cData = [cData C];
    end
end