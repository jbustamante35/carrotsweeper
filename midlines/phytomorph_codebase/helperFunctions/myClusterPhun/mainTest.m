%%
clear cT
clear class clusterTree clusterNode clusterForest arborist
%%


% creaete extract function to give to the arborist
filter = fspecial('gaussian',[21 21],5);
extractFunc = @(X)rgbAndhsvExtract(X,filter);
% create tree growth rules for the arborist
suggestNumberOfClusterFunction = @(data,level,treePara)treePara(1)*(level<=treePara(2)) + 1*(level>treePara(2));
% create cluster parameter generating function for arborist
clusterFunctionFunctionGenerator = @(X,K)fitgmdist(X,K,'Options',options,'Start','plus','RegularizationValue',0.0001);


maxBranch = 4;
maxDepth = 2;



jA = arborist(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc,maxDepth,maxBranch);

forest = jA.plantTrees(HISTO_TOT,3);
%%
k = forest.treeSet{1}.clusterImage(FileList{1}{5});
%%
k = forest.clusterImage(FileList{1}{5});
%%

jM = johnMuir();
metricFunc = @(X)blobMeasure_redCaps(X);
jM.attachMetric(metricFunc);

iT = jM.measureTree(cT,FileList{5}(2));
%%
iF = jM.measureForest(forest,FileList{5}(2),2);
%%


options = statset('Display','iter');
clusterFunctionFunctionGenerator = @(X,K)fitgmdist(X,K,'Options',options,'Start','plus','RegularizationValue',0.0001);
cT = clusterTree(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc);
cT.build(HISTO_TOT);





k = cT.clusterImage(FileList{1}{5},extractFunc,filter);
imshow(k,[]);


[m] = blobMeasure_redCaps(k==4);


close all
UQ = unique(k);
for u = 1:numel(UQ)
    imshow(k==UQ(u),[]);
    drawnow
    waitforbuttonpress
end