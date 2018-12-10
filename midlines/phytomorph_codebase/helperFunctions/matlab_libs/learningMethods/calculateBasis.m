function [basisInformation] = calculateBasis(selectFunction,classifyFeatures,transformFeatures,domainFeatureMap,codomainFeatureMap)
    % select points from domain map
    [indexSets] = featureMapBank.inverseSelect(selectFunction,classifyFeatures,transformFeatures,domainFeatureMap);
    % intersect the domain index with the codomain index
    [indexSets] = featureMapBank.intersectIndex(indexSets,{codomainFeatureMap});
    % load features from the index intersection
    [data] = featureMapBank.loadFeatureMapsAtAcrossIndexSets({codomainFeatureMap},indexSets);
    % decompose
    [S C U E L ERR LAM] = PCA_FIT_FULL_T(data,size(data,1));
    % store and return
    basisInformation.E = E;
    basisInformation.U = U;
end