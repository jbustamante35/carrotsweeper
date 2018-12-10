FilePath = '/mnt/tetra/nate/overHead/';
FileList = {};
FileExt = {'jpg','jpeg'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% remove dark images
for s = 1:numel(SET)
    value = [];
    for e = 1:numel(SET{s})
        I = imread(SET{s}{e});
        G = rgb2gray(I);
        value(e) = mean(G(:));
    end
    SET{s}(value < 50) = [];
    s
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create arborist for building trees
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
SAM = jA.sampleElements(SET,[2 30 40000 .5]);
%% plant forest
forest = jA.plantTrees(SAM,5);
%%
k = forest.clusterImage(SET{1}{end});
%%
for e = 1:5
    oI = imread(SET{e}{end});
    [k p] = forest.clusterImage(oI);
    for i = 1:size(k,3)
        rgb = label2rgb(k(:,:,i));
        sig = cat(2,oI,rgb);
        imshow(sig,[]);
        title(num2str(i))
        waitforbuttonpress
    end
end
%%
overHeadFunc = @(X)extractPhenotypesFromOverheadCamera_ver2(X,forest.getTree(1),'./output/');
%overHeadFunc(SET{1})
%%
pF = partialFunction(overHeadFunc,'overHeadFuncAPP');
pF.publish();
% run local =- cornPopper(FileList{91})
%%
close all
for e = 1%:size(k,3)
    imshow(k(:,:,e),[]);
    drawnow
    waitforbuttonpress
end
%% cluster
K = 5;
SAM = 5;
options = statset('Display','iter','MaxIter',500);
%levelD = fitgmdist(kC(:,[1]),K,'Start','plus','Options',options);
levelD = fitgmdist(CS(1:SAM:end,:),K,'Start','plus','Options',options,'RegularizationValue',.000001);
%% view cluster level 1
s = 1;
e = 1;
oI = double(imread(SET{s}{e}));
osz = size(oI);
I = permute(oI,[3 1 2]);
sz = size(I);
I = reshape(I,[sz(1) prod(sz(2:3))]);
idx = cluster(levelD,double(I'));
idx = reshape(idx,osz(1:2));
close all
imshow(idx,[])
%% sub cluster
K =2;
tempyIDX = cluster(levelD,CS);
SAM = 5;
options = statset('Display','iter');
levelDD = fitgmdist(CS(tempyIDX==2,:),K,'Start','plus','Options',options);
%% sub cluster
K = 3;
tempyIDX = cluster(levelD,CS);
SAM = 5;
options = statset('Display','iter');
levelDD = fitgmdist(CS(tempyIDX==1,:),K,'Start','plus','Options',options);

%%
oPath = './cullenLocal/testForCullen_againX/';
for s = 1:numel(SET)
    tmpPath = strrep(oPath,'X',num2str(s));
    phMeasure = extractPhenotypesFromOverheadCamera_ver2(SET{s},forest.getTree(1),tmpPath);
end
%%
overHeadFunc= @(X)extractPhenotypesFromOverheadCamera_ver2(X,forest.getTree(1),'./output/');
pF = partialFunction(overHeadFunc,'overHeadFuncAPP');
pF.publish();
% run local =- cornPopper(FileList{91})
%%



CS = [CS;I'];