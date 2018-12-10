function [ret] = pbPCA_loop(I,points,domainPara)
    % generate domain
    tD = phytoDomain(domainPara);
    tD.generateDomain();
    % loop over points in list [dims, numberPoints]
    for pt = 1:size(points,2)
        % get current point
        curPoint = points(:,pt);
        % make affine
        tA = eye(3);
        tA(1:2,3) = curPoint;
        % point balanced pca
        Tout = pbPCA(I,tA,tD.d,@lt);
        % make affine
        Tout(1:2,3) = curPoint;
        % make return
        ret{pt} = phytoAffine(Tout);
    end
end