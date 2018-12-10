FilePath = '/mnt/tetra/nate/mwald/';
FileList = {};
FileExt = {'jpg'};
FileList = sdig(FilePath,FileList,FileExt,1);
%% find the straights
kidx = [];
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e}{1},'SG'))
        kidx = [kidx e];
    end
end
sFileList = FileList(kidx);
%% load the straight growth
T = linspace(-pi,pi,100)*180/pi;
for e = 1:numel(sFileList)
    tmpI = imread(sFileList{e}{1});
    tmpI = imresize(tmpI,.25);
    for r = 1:numel(T)
        tmpR = imrotate(tmpI,T(r),'crop');
        imshow(tmpR,[]);
        drawnow
    end
end
%% load first image from each set
I = [];
for e = 1:numel(FileList)
    I(:,:,:,e) = imread(FileList{e}{1});
    imshow(I(:,:,:,e)/255,[])
    drawnow
end
%% 
%% try to learn if rotates
FilePath = '/mnt/tetra/nate/mwald/';
FileList = {};
FileExt = {'jpg'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
S = [];
parfor e = 1:numel(FileList)
    tmp = imread(FileList{e});
    tmp = imresize(tmp,.25);
    tmp = rgb2gray(tmp);
    S(:,:,e) = tmp;
    if ~isempty(strfind(FileList{e},'SG'))
        Y{e} = 'straight';
    else
        Y{e} = 'non-straight';
    end
    fprintf(['done loading:' num2str(e) ':' num2str(numel(FileList)) '\n']);
end
%% loop over stack
close all
for e = 1:size(S,3)
    imshow(S(:,:,e),[]);
    drawnow
end
%% crop ove stack
[J,BOX] = imcrop(S(:,:,1)/255);
cS = [];
for e = 1:size(S,3)
    cS(:,:,e) = imcrop(S(:,:,e),BOX);
    fprintf(['done cropping:' num2str(e) ':' num2str(size(S,3)) '\n']);
end
%% look at cropped stack
for e = 1:size(cS,3)
    imshow(cS(:,:,e),[]);
    drawnow
end
%% make inputs the right size
sz = size(cS);
cS = reshape(cS,[sz(1) sz(2) 1 sz(3)]);
%% train network
layers = [imageInputLayer([size(cS,1) size(cS,2) 1]);
          convolution2dLayer([21 21],3);
          reluLayer();
          maxPooling2dLayer([10 10],'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
options = trainingOptions('sgdm','MaxEpochs',2,'InitialLearnRate',0.0001,'ExecutionEnvironment','parallel');
rotNet = trainNetwork(cS,categorical(Y'),layers,options);
%%
Ypre = classify(rotNet,cS);



