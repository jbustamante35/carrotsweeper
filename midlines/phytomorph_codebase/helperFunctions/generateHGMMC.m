function [hgmm] = generateHGMMC(H,clusterLevels)
    options = statset('Display','iter');
    hgmm(1).gmm = fitgmdist(H,clusterLevels(1),'Options',options,'Start','plus','RegularizationValue',0.0001);
    
    
    
    idx = cluster(hgmm(1).gmm,H);
    for l = 1:clusterLevels(1)
        if numel(clusterLevels) ~= 1
            [hgmm(1).nextLevel(l)] = generateHGMMC(H(idx==l,:),clusterLevels(2:end));
        end
    end
    
end