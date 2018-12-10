function [S C U E ERR] = clusterPCA(data,kidx,percent)
    S = [];
    ERR = [];
    % find the data points in the cluster
    UQ = unique(kidx);
    % find the data points in the cluster
    for u = 1:numel(UQ)
        % find the data points in the cluster
        fidx = find(kidx == UQ(u));
        % get the sub data
        sdata = data(fidx,:);
        % perform all dims latent 
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(sdata,size(sdata,2));
        
        % find dims for percent explained
        dims = find(tL < percent);
        dims = max(dims);
        
        % find latent model for percent explained
        [fS C{u} U{u} E{u} fL fERR fLAM] = PCA_FIT_FULL(sdata,dims);
        % stack SIM and error
        S = [S;fS];
        ERR = [ERR;fERR];
    end
end