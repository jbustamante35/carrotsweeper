FilePath = '/home/nate/Downloads/scottbrainard/';
FileList = {};
FileExt = {'NEF'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% sample NEFS
% init the arborist with rules for building trees
options = statset('Display','iter');
% creaete extract function to give to the arborist
filter = fspecial('gaussian',[21 21],5);
extractFunc = @(X)rgbAndhsvExtract(X,filter);
% create tree growth rules for the arborist
suggestNumberOfClusterFunction = @(data,level,treePara)treePara(1)*(level<=treePara(2)) + 1*(level>treePara(2));
% create cluster parameter generating function for arborist
clusterFunctionFunctionGenerator = @(X,K)fitgmdist(X,K,'Options',options,'Start','plus','RegularizationValue',0.0001);
% spec the tree parameters
maxBranch = 4;
maxDepth = 3;
% build the arborist
jA = arborist(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc,maxDepth,maxBranch);
%% sample data with johnny
SAM = jA.sampleElements(FileList,[40 1 30000 .25]);
%%
% plant the trees
forest = jA.plantTrees(SAM,10);
%%
close all
for i = 1:30
    oI = double(imread(FileList{i}));
    k = forest.clusterImage(oI);
    for e = 1:size(k,3)
        rgb = label2rgb(k(:,:,e));
        imshow(cat(2,oI/255,double(rgb)/255),[]);
        drawnow
        waitforbuttonpress
    end
end