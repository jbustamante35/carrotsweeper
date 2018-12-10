function [ret] = pbPCA_file_loop(I,objSet,tD,timeIndex)
    % loop over objects in objectSet
    for pt = 1:numel(objSet)
        
        % get current object in the set
        curObj = objSet{pt};
        
        % the the time index from the point
        curPoint = curObj(:,timeIndex);
        
        % the the time index from the point
        tA = eye(3);
        tA(1:2,3) = curPoint(1:2);
        
        % point balanced pca
        Tout = pbPCA(I,tA,tD.d,@lt);
        
        % make affine
        Tout(1:2,3) = curPoint(1:2);
        % make return
        ret{pt} = phytoAffine(Tout);
    end
    
end


%{




% generate domain
    tD = phytoDomain(domainPara);
    tD.generateDomain();
    
    % loop over points in list [dims, numberPoints]
    for pt = 1:numel(geoSet)
        
        % get current object
        curObj = geoSet{pt};
        % get point representation of object
        
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
%}