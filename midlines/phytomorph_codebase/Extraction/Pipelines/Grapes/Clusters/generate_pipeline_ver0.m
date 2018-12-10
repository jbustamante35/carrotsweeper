FilePath = '/home/nate/Downloads/grapes/other/GE1025 2017 Color Corrected TIFF/';
FileList = {};
FileExt = {'tif'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
re = .5;
S = {};
cnt = 1;
for e = 1:5:numel(FileList)
    I = imread(FileList{e});
    I = imresize(I,re);
    I = permute(I,[3 1 2]);
    sz = size(I);
    I = reshape(I,[sz(1) prod(sz(2:3))]);
    S{cnt} = I;
    cnt = cnt + 1;
    fprintf(['Done with :' num2str(cnt) '\n']);
end
%% stack from cell
CD = [];
for e = 1:numel(S)
    CD = [CD;S{e}'];
end
%%
CD = CD(1:5:end,:);
%%
CD = double(CD);
%%
% init the arborist with rules for building trees
options = statset('Display','iter');
% creaete extract function to give to the arborist
filter = fspecial('gaussian',[21 21],5);
extractFunc = @(X)rgbAndhsvExtract(X,filter);
% create tree growth rules for the arborist
suggestNumberOfClusterFunction = @(data,level,treePara)treePara(1)*(level<=treePara(2)) + 1*(level>treePara(2));
% create cluster parameter generating function for arborist
clusterFunctionFunctionGenerator = @(X,K)fitgmdist(X,K,'Options',options,'Start','plus','RegularizationValue',0.0001);
% feature selection function
idxSelectorFunction = @(X,L)logical([1 1 1 1 1]);
% spec the tree parameters
maxBranch = 3;
maxDepth = 3;
% build the arborist
jA = arborist(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc,maxDepth,maxBranch,idxSelectorFunction);
%% sample data with johnny
SAM = jA.sampleElements(FileList,[400 1 10000 .5]);
%% plant the trees
forest = jA.plantTrees(SAM,10);
%% tone down filter
filter = fspecial('gaussian',[21 21],1);
extractFunc = @(X)rgbAndhsvExtract(X,filter);
for e = 1:forest.numberOfTrees
    forest.treeSet{e}.extractFunc = extractFunc;
end
%%
close all
for i = 1:30
    oI = double(imread(FileList{i}));
    [k] = forest.clusterImage(FileList{i});
    for e = 1:size(k,3)
        rgb = label2rgb(k(:,:,e));
        imshow(cat(2,oI/255,double(rgb)/255),[]);
        drawnow
        waitforbuttonpress
    end
end
%%
%% test out john
jM = johnMuir();
metricFunc = @(X)grapeCluster_MaskMeasure(X);
jM.attachMetric(metricFunc);
overlayFunc = @(I,M,fileName)generateGrapeOverLay(I,M,fileName);
jM.attachOverlayFunction(overlayFunc,'./output/grapeTrees2/');
%%
close all
iF = jM.measureForest(forest,FileList(1:10),2);
%% isolate trees
close all
iF = jM.measureForest(forest,FileList(1:50),3,[3 7]);
%%
for i = 1:10
    masterMeasure_ver0(FileList{i},forest,'./output/');
end


