function [] = manifold_learn(X,Y,para)
    % decompose via linear
    [S C U E L ERR LAM] = PCA_FIT_FULL(X,para.comp);
    % split into groups
    kidx = kmeans(C,para.groupN);
    % loop over each group
    UQ = unique(kidx);
    for u = 1:numel(UQ)
        subX = X(kidx==UQ(u),:);
        [S C U E L ERR LAM] = PCA_FIT_FULL(X,para.comp);
        subY = Y(kidx==UQ(u),:);
        pVec = subY/subX;
    end
end