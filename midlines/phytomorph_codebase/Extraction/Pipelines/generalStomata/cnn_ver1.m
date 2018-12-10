%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIZE PIPELINE - start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Color images - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gather and color images for maize
%% scan for maize NMS files - RILS
FilePath = '/mnt/tetra/nate/RILs/';
iFileList = {};
FileExt = {'nms'};
iFileList = gdig(FilePath,iFileList,FileExt,1);
%% gather clicks for training data stomata Centers
for e = 1:30
    I = imread(iFileList{e});
    [maizeStomataCenter_row{e} maizeStomataCenter_column{e} v{e}] = impixel(I);
    e
end
%% gather clicks for training data - stomata area
BOX_size = [80 40];
for e = 1:30
    I = imread(iFileList{e});
    for p = 1:numel(maizeStomataCenter_row{e})
        BOX = [maizeStomataCenter_row{e}(p) - BOX_size(1)/2 maizeStomataCenter_column{e}(p) - BOX_size(2)/2 ...
            BOX_size];
        subI = imcrop(I,BOX);
        imshow(subI,[]);
        drawnow
    end
end
%% save clicks
%save('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/RIL_maize_stomata_centers.mat','maizeStomataCenter_row','maizeStomataCenter_column');
%% load clicks
%load('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/RIL_maize_stomata_centers.mat','maizeStomataCenter_row','maizeStomataCenter_column');
%% color images from clicks
I = imread(iFileList{1});
border = 40;
close all
Y = [];
for e = 1:numel(maizeStomataCenter_row)
    tmp = zeros(size(I));
    %tmpI = imread(iFileList{e});
    for p = 1:size(maizeStomataCenter_row{e},1)
         tmp(maizeStomataCenter_column{e}(p),maizeStomataCenter_row{e}(p)) = 1;
    end
    tmp = imdilate(tmp,strel('disk',11,0));
    for rot = 1:4
        tmp(1:border,:) = [];
        %tmpI(1:border,:) = [];
        tmp = imrotate(tmp,90);
        tmpI = imrotate(tmpI,90);
    end
    tmp = padarray(tmp,[1 1],0,'pre');
    tmpI = padarray(tmpI,[1 1],0,'pre');
    size(tmp)
    
    %out = flattenMaskOverlay(tmpI,logical(tmp));
    %imshow(out,[]);
    Y = [Y;tmp(:)];
    
    
    
    drawnow
end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Color images - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract images - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% extract bug eye
for e = 1:30
    [fft{e}] = extractSingleBugEye_v2(iFileList{e},15,[1 40]);
end

%% stack data
imgTOT = 433*433;
RAD = size(fft{1}.f,3);
TH = size(fft{1}.f,4);
X = zeros([RAD TH 3 imgTOT*numel(fft)]);
str = 1;
for e = 1:numel(fft)
    
    fprintf(['start stacking:' num2str(e) ':' num2str(numel(fft)) '\n']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp = permute(fft{e}.f,[3 4 1 2]);
    tsz = size(tmp);
    tmp = reshape(tmp,[tsz(1:2) 433*433]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpA = fft{e}.A;
    tmpA = diff(cat(4,tmpA,tmpA(:,:,:,end)),1,4);
    tmpA = abs(tmpA)/(2*pi);
    tmpA = .5*(((tmpA.^2)+1).^.5 + (((1-tmpA).^2)+1).^.5);
    tmpA = permute(tmpA,[3 4 1 2]);
    tmpA = reshape(tmpA,[tsz(1:2) 433*433]);
   
    
    
    stp = str + size(tmpA,3) - 1;
    X(:,:,:,str:stp) = cat(3,reshape(tmp,[tsz(1:2) 1 prod(tsz(3:4))]),...
                             reshape(tmpA,[tsz(1:2) 1 prod(tsz(3:4))]),...
                             zeros([tsz(1:2) 1 prod(tsz(3:4))]));
    
    str =  stp + 1;
    fprintf(['end stacking:' num2str(e) ':' num2str(numel(fft)) '\n'])
end
%% convert to single
X = single(X);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract images - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Traing CNN - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define CNN
layers = [imageInputLayer([40 20 3])
          convolution2dLayer([3,3],10)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          convolution2dLayer([3,3],10)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
      
layers = [imageInputLayer([40 20 3])
          convolution2dLayer([7,2],10)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
      
layers = [imageInputLayer([size(X,1) size(X,2) 3])
          convolution2dLayer([7,2],20)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          convolution2dLayer([3,3],10)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
      
NTR = round(1*numel(Y));
options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3);
%% train CNN
convnet_Maize2 = trainNetwork(X(:,:,:,1:NTR),categorical(Y(1:NTR)),layers,options);
%save('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/maize_CNN.mat','convnet_Maize');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Traing CNN - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Run local and condor - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather data from CyVerse
% this is the RILs only path
%dataPath = '/iplant/home/leakey_cyverse/maizeData/stomataTopoData/Hybrid';
%dataPath = '/iplant/home/leakey_cyverse/maizeData/stomataTopoData/Inbred';
dataPath = '/iplant/home/leakey_cyverse/maizeData/stomataTopoData/RILs';
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
FileList = {};
FileExt = {'nms'};
 for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
 end
%% search for file
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'601008 leaf2-2'))
        e
    end 
end
%% issue tickets over the FileList
[FileList] = issueBulkTicket(FileList);
%% issue ticket over the return folder
%remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_maizeHybrid_ver0/'];
remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_maizeRIL_ver1/'];
[remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numel(FileList),'write');
%% run test local
stomata_cnnVersion(FileList{1000},convnet_Maize2,15,[1 40],'','');
%% lauch on condor
func = cFlow('stomata_cnnVersion');
func.setMCRversion('v920');
func.setMemory('8000');
for e = 1:numel(FileList)
    func(FileList{e},convnet_Maize,15,[1 40],'./output/',remoteOutputLocation);
    fprintf(['done rendering job:' num2str(e) ':' num2str(numel(FileList)) '\n'])
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
    
func.submitDag(auth,150,150);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Run local and condor - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIZE PIPELINE - end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORGHUM PIPELINE - start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Color images - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather data from CyVerse
% this is the RILs only path
dataPath = '/iplant/home/leakey_cyverse/sorghumData/stomataTopoData/Accessions_2016';  
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
sorghum_FileList = {};
FileExt = {'nms'};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        sorghum_FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%% gather clicks for training data stomata Centers
for e = 1:30
    I = imread(sorghum_FileList{e});
    [sorghumStomataCenter_row{e} sorghumStomataCenter_column{e} v{e}] = impixel(I);
    e
end

%{
%% gather clicks for training data - stomata area
BOX_size = [80 40];
for e = 1:30
    I = imread(iFileList{e});
    for p = 1:numel(maizeStomataCenter_row{e})
        BOX = [maizeStomataCenter_row{e}(p) - BOX_size(1)/2 maizeStomataCenter_column{e}(p) - BOX_size(2)/2 ...
            BOX_size];
        subI = imcrop(I,BOX);
        imshow(subI,[]);
        drawnow
    end
end
%}
%% save clicks
%save('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/RIL_sorghum_stomata_centers.mat','sorghumStomataCenter_row','sorghumStomataCenter_column');
%% load clicks
load('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/RIL_sorghum_stomata_centers.mat','sorghumStomataCenter_row','sorghumStomataCenter_column');
%% color images from clicks
I = imread(sorghum_FileList{1});
border = 20;
close all
Y = [];
for e = 1:numel(sorghumStomataCenter_column)
    I = imread(sorghum_FileList{e});
    
    tmp = zeros(size(I));
    for p = 1:size(sorghumStomataCenter_column{e},1)
         tmp(sorghumStomataCenter_column{e}(p),sorghumStomataCenter_row{e}(p)) = 1;
    end
    
    
    tmp2 = imdilate(tmp,strel('disk',9));
    tmp = imdilate(tmp,strel('disk',5,0));
    tmp = tmp2 + tmp;
    tmp = tmp2;
    
    for rot = 1:4
        tmp(1:border,:) = [];
        I(1:border,:) = [];
        tmp = imrotate(tmp,90);
        I = imrotate(I,90);
    end
    tmp = padarray(tmp,[1 1],0,'pre');
    
    I = padarray(I,[1 1],0,'pre');
    
    
    out = flattenMaskOverlay(bindVec(I),logical(tmp));
    
    
    imshow(out,[]);
    drawnow
    Y = [Y;tmp(:)];
    
    
    
    drawnow
end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Color images - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract images - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% extract bug eye
% radius was 40 but bad job me changed to 20
for e = 1:30
    fprintf(['start extracting:' num2str(e) '\n']);tic
    [fft{e}] = extractSingleBugEye_v2(sorghum_FileList{e},15,[3 20],[]);
    fprintf(['done extracting:' num2str(e) ':' num2str(toc) '\n'])
end
%% extract im2col style
pSZ = 33;
I = imread(sorghum_FileList{1});
I = im2col(I,[pSZ pSZ],'sliding');
X = zeros(pSZ,pSZ,1,size(I,2)*30);
Z = size(I,2);
str = 1;
diskSample = 1;

if diskSample
    [R,T] = ndgrid(linspace(0,(pSZ-1)/2,(pSZ-1)/2),linspace(-pi,pi,round(2*pi*(pSZ-1)/2)));
    X1 = R.*cos(T) + (pSZ-1)/2;
    X2 = R.*sin(T) + (pSZ-1)/2;
    X = zeros([size(X1) size(X,3) size(X,4)]);
end


for e = 1:30
    fprintf(['start extracting:' num2str(e) '\n']);tic
    I = imread(sorghum_FileList{e});
    I = im2col(I,[pSZ pSZ],'sliding');
    I = reshape(I,[pSZ pSZ 1 size(I,2)]);
    
    if diskSample
        nI = zeros([size(X1) size(I,3) size(I,4)]);
        parfor s = 1:size(nI,4)
            nI(:,:,:,s) = ba_interp2(I(:,:,:,s),X1,X2);
        end
        I = nI;
    end
    
    
    stp = str + Z - 1;
    X(:,:,:,str:stp) = I;
    str = stp + 1;
    fprintf(['done extracting:' num2str(e) ':' num2str(toc) '\n'])
end
%% resample via disk
[R,T] = ndgrid(linspace(0,pSZ,pSZ),linspace(pi,pi,round(2*pi*pSZ)));
X1 = R.*cos(T) + (pSZ-1)/2;
X2 = R.*sin(T) + (pSZ-1)/2;
nX = zeros([size(X1) size(X,3) size(X,4)]);
parfor e = 1:size(X,4)
    nX(:,:,:,e) = ba_interp2(X(:,:,:,e),X1,X2);
    e
end
%% stack data
imgTOT = 473*473;
RAD = size(fft{1}.f,3);
TH = size(fft{1}.f,4);
X = zeros([RAD TH 3 imgTOT*numel(fft)]);
str = 1;
for e = 1:numel(fft)
    
    fprintf(['start stacking:' num2str(e) ':' num2str(numel(fft)) '\n']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp = permute(fft{e}.f,[3 4 1 2]);
    tsz = size(tmp);
    tmp = reshape(tmp,[tsz(1:2) prod(tsz(3:4))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpA = fft{e}.A;
    tmpA = diff(cat(4,tmpA,tmpA(:,:,:,end)),1,4);
    tmpA = abs(tmpA)/(2*pi);
    tmpA = .5*(((tmpA.^2)+1).^.5 + (((1-tmpA).^2)+1).^.5);
    tmpA = permute(tmpA,[3 4 1 2]);
    tmpA = reshape(tmpA,[tsz(1:2) prod(tsz(3:4))]);
   
    
    
    stp = str + size(tmpA,3) - 1;
    X(:,:,:,str:stp) = cat(3,reshape(tmp,[tsz(1:2) 1 prod(tsz(3:4))]),...
                             reshape(tmpA,[tsz(1:2) 1 prod(tsz(3:4))]),...
                             zeros([tsz(1:2) 1 prod(tsz(3:4))]));
    
    str =  stp + 1;
    fprintf(['end stacking:' num2str(e) ':' num2str(numel(fft)) '\n'])
end
%% convert to single
X = single(X);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract images - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Traing CNN - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define CNN
layers = [imageInputLayer([size(X,1) size(X,2) 3])
          convolution2dLayer([1,5],15)
          reluLayer
          %maxPooling2dLayer([2 2],'Stride',2)
          convolution2dLayer([7,1],10)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];

layers = [imageInputLayer([size(X,1) size(X,2) 3])
          convolution2dLayer([1,5],5) % change from 10 to 5
          reluLayer
          maxPooling2dLayer([1 2],'Stride',1)
          convolution2dLayer([1,5],5) % change from 10 to 5
          reluLayer
          convolution2dLayer([7,2],5)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];

      
NTR = round(1*numel(Y));
options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3,'Minibatch',128*4,'CheckpointPath','/mnt/snapper/nate/CP/','Plots','training-progress');
final_i_hope = trainNetwork(X,categorical(Y),layers,options);
%% simple test
CC = predict(final_i_hope,X(:,:,:,1:473*473));
TT = Y(1:473*473);
I = imread(sorghum_FileList{1});
%% apply FFT to X
fX = X;
for e = 1:size(X,4)
    fX(:,:,:,e) = fft2(X(:,:,:,e));
    e
end
%% diff rot invar
SKIP = 2;
layers = [imageInputLayer([size(X,1) size(X,2) 1],'Normalization','None')
          convolution2dLayer([1,45],15,'padding','same')
          convolution2dLayer([9,1],4,'padding','same')
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          convolution2dLayer([2,2],3,'Padding','same')
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];

      
layers = [imageInputLayer([size(X,1) size(X,2) 1],'Normalization','None')
          convolution2dLayer([1,5],5,'padding','same')
          reluLayer
          maxPooling2dLayer([1,2],'Stride',2)
          convolution2dLayer([1,5],5,'padding','same')
          convolution2dLayer([5,1],3,'padding','same')
          reluLayer
          convolution2dLayer([2,2],3,'Padding','same')
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
%imageAugmenter = imageDataAugmenter('RandRotation',[-10 10],'RandXScale',[.6 1.2],'RandYScale',[.6 1.2]);
%imageAugmenter = imageDataAugmenter('RandRotation',[-90 90]);
%imageSize = [33 33 1];
imageSize = [size(X,1) size(X,2)];
NTR = 4;
TR = 25;
VAL = (30 - TR);
strTR = 1;
dSZ = 472;
stpTR = dSZ*dSZ*TR;
strVAL = stpTR + 1;
stpVAL = numel(Y);

%(:,:,:,3*230400:end)
%VALData = {X(:,:,:,strVAL:NTR:stpVAL),categorical(Y(strVAL:NTR:stpVAL))};
%datasource = augmentedImageSource(imageSize,X,categorical(Y),'DataAugmentation',imageAugmenter);
datasource = augmentedImageSource(imageSize,X(:,:,:,3*230400:end),categorical(Y(3*230400:end)));
%options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3,'CheckpointPath','/mnt/snapper/nate/CP/','Plots','training-progress','ValidationFrequency',200,'ValidationData',VALData);
options = trainingOptions('sgdm','InitialLearnRate',.001,'Minibatch',128*8,'ExecutionEnvironment','parallel','MaxEpochs',3,'CheckpointPath','/mnt/snapper/nate/CP/','Plots','training-progress');
%convnet_Sorghum6 = trainNetwork(datasource,convnet_Sorghum6.Layers,options);
convnet_Sorghum8_3s = trainNetwork(datasource,layers,options);
%convnet_Sorghum6 = trainNetwork(datasource,layers,options);
%% expand and train
[newLayers] = expandNet(convnet_Sorghum8_3s);
convnet_Sorghum9_3s = trainNetwork(datasource,newLayers,options);
%% view #0
close all
fidxSEL = find(Y==1);
Features = activations(convnet_Sorghum8_2,X(:,:,1,fidxSEL(end-100)),'conv_1','OutputAs','channels');
sz = size(Features);
Features = reshape(Features,[sz(1) sz(2) 1 sz(3)]);
montage(mat2gray(Features),'Size',[4 5]);
%% vew #1
IT = imread(sorghum_FileList{3521});
[c r V] = impixel(IT,[]);
%% view #2
close all
subI = IT(r-16:r+16,c-16:c+16);
features = activations(convnet_Sorghum8_2s,subI,'conv_1');
features = reshape(features,[33 33 5]);
[maxValue,maxValueIndex] = max(max(max(abs(features))));
imshow(cat(2,subI,features(:,:,maxValueIndex)),[]);
featuresD = activations(convnet_Sorghum8_2s,subI,'conv_2','OutputAs','channels');
[maxValue5,maxValueIndex5] = max(max(max(abs(featuresD))));
act5chMax = featuresD(:,:,maxValueIndex5);
figure;
imshow(cat(1,subI,imresize(abs(act5chMax),size(subI))),[]);
%%
[newLayers] = expandNet(convnet_Sorghum8);
convnet_Sorghum9 = trainNetwork(datasource,newLayers,options);
%%
layers = [imageInputLayer([size(X,1) size(X,2) 1])
          convolution2dLayer([9,7],9,'Padding',round((7-1)/2))
          reluLayer
          maxPooling2dLayer([2 2],'Stride',1)
          convolution2dLayer([3,3],3,'Padding',round((3-1)/2))
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
%% test multi gpu
func = cFlow('trainType1');
func.setMCRversion('v930');
func.setGPU(4);
func.setMemory('35000');
layers = [imageInputLayer([size(X,1) size(X,2) 1])
      convolution2dLayer([7,7],9,'Padding',round((7-1)/2))
      reluLayer
      maxPooling2dLayer([2 2],'Stride',1)
      convolution2dLayer([4,4],3,'Padding',round((7-1)/2))
      fullyConnectedLayer(2)
      softmaxLayer
      classificationLayer()];
  
net = func(X,Y,layers,1);


auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%%
cTR = cvpartition(Y,'KFold',2);


optimVars = [
    optimizableVariable('NetworkDepth',[1 3],'Type','integer')
    optimizableVariable('filterSize',[3 9],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
maxE = 1;
toSlow = 1;
exeEnvironment = {'cpu',false};

%{
func = cFlow('trainType1');
func.setMCRversion('v930');
func.setGPU(2);
func.setMemory('20000');

func = cFlow('hyperPdeploy');
%}

[BayesObject] = hyperPdeploy(optimVars,X(:,:,:,find(cTR.test(1))),...
                            categorical(Y(find(cTR.test(1)))),...
                            X(:,:,:,find(cTR.test(2))),...
                            categorical(Y(find(cTR.test(2)))),...
                            exeEnvironment,maxE,toSlow,2);
%%
[valError,cons,NETfileName] = makeObjFcn(me,X,categorical(Y),[],[],0,'training-progress',1,'parallel',2);
                        %%
    BayesObject = bayesopt(@(OPS)ObjFcn(OPS,X(:,:,:,find(cTR.test(1))),categorical(Y(find(cTR.test(1)))),X(:,:,:,find(cTR.test(2))),categorical(Y(find(cTR.test(2))))),...
        optimVars,...
        'MaxObj',90,...
        'MaxTime',2*60*60,...
        'IsObjectiveDeterministic',false,...
        'UseParallel',true);
%%
%save('/home/nate/b.mat','BayesObject');
load('/home/nate/b.mat','opt');
[options,Layers] = autoBuildNetwork(opt,'training-progress',8,numel(unique(Y)));
imageAugmenter = imageDataAugmenter('RandRotation',[-90 90],'RandXScale',[.9 1.1],'RandYScale',[.9 1.1]);
datasource = augmentedImageSource([41 41 1],X,categorical(Y),'DataAugmentation',imageAugmenter);
trainedNet = trainNetwork(datasource,Layers,options);
%% 
load('/mnt/tetra/nate/CNN/convnet_checkpoint__12115__2017_10_26__21_53_45.mat');
%% train CNN
convnet_Sorghum5 = trainNetwork(X(:,:,:,1:NTR),categorical(Y(1:NTR)),layers,options);
%save('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/sorghum_CNN4.mat','convnet_Sorghum4');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Traing CNN - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Run local and condor - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather data from CyVerse
% this is the RILs only path
dataPath = '/iplant/home/leakey_cyverse/sorghumData/stomataTopoData/Accessions_2016'; 
%dataPath = '/iplant/home/leakey_cyverse/CharlesPignonData';
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
sorghum_FileList = {};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        sorghum_FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%% search for file
for e = 1:numel(sorghum_FileList)
    if ~isempty(strfind(sorghum_FileList{e},'EF0175_2_1'))
        e
    end 
end
%% issue tickets over the FileList
[sorghum_FileList] = issueBulkTicket(sorghum_FileList);
%% issue ticket over the return folder
remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_sorghum2016_verFinal11/'];
%remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_forCharles/'];
[remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numel(sorghum_FileList),'write');
%% run test local
figure
stomata_cnnVersion(sorghum_FileList{3521},convnet_Sorghum5,15,[3 20],'','');
%% run test local
figure
%stomata_cnnVersion3(sorghum_FileList{3521},convnet_Sorghum8_2s,15,[3 20],'','');
%stomata_cnnVersion3(sorghum_FileList{3159},convnet_Sorghum8_2s,15,[3 20],'','');
%stomata_cnnVersion4(sorghum_FileList{1752},final_i_hope,15,[3 20],'','',[]);
stomata_cnnVersion(sorghum_FileList{1752},final_i_hope,15,[3 20],'','');
%stomata_cnnVersion3(sorghum_FileList{352},convnet_Sorghum9_2,15,[3 20],'','');
%stomata_cnnVersion3(sorghum_FileList{5643},convnet_Sorghum6,15,[3 20],'','');

%% launch on condor
func = cFlow('stomata_cnnVersion3');
func.setMCRversion('v930');
func.setMemory('4000');

toRun = numel(sorghum_FileList);
%toRun = 5;
for e = 1:toRun
    func(sorghum_FileList{e},convnet_Sorghum9_2s,15,[3 40],'./output/',remoteOutputLocation);
    fprintf(['done rendering job:' num2str(e) ':' num2str(toRun) '\n'])
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
    
func.submitDag(auth,250,250);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Run local and condor - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORGHUM PIPELINE - end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETARIA PIPELINE - start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Color images - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather data from CyVerse
% this is the RILs only path
dataPath = '/iplant/home/leakey_cyverse/setariaData/stomataTopoData%';  
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
setaria_FileList = {};
FileExt = {'nms'};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        setaria_FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%% gather clicks for training data stomata Centers
for e = 1:30
    I = imread(setaria_FileList{e});
    [setariaStomataCenter_row{e} setariaStomataCenter_column{e} v{e}] = impixel(I);
    e
end

%{
%% gather clicks for training data - stomata area
BOX_size = [80 40];
for e = 1:30
    I = imread(iFileList{e});
    for p = 1:numel(maizeStomataCenter_row{e})
        BOX = [maizeStomataCenter_row{e}(p) - BOX_size(1)/2 maizeStomataCenter_column{e}(p) - BOX_size(2)/2 ...
            BOX_size];
        subI = imcrop(I,BOX);
        imshow(subI,[]);
        drawnow
    end
end
%}
%% save clicks
%save('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/setaria_stomata_centers.mat','setariaStomataCenter_row','setariaStomataCenter_column');
%% load clicks
%load('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/setaria_stomata_centers.mat','setariaStomataCenter_row','setariaStomataCenter_column');
%% color images from clicks
I = imread(setaria_FileList{1});
border = 40;
close all
Y = [];
for e = 1:numel(setariaStomataCenter_column)
    I = imread(setaria_FileList{e});
    
    tmp = zeros(size(I));
    for p = 1:size(setariaStomataCenter_column{e},1)
         tmp(setariaStomataCenter_column{e}(p),setariaStomataCenter_row{e}(p)) = 1;
    end
    tmp = imdilate(tmp,strel('disk',11,0));
    for rot = 1:4
        tmp(1:border,:) = [];
        I(1:border,:) = [];
        tmp = imrotate(tmp,90);
        I = imrotate(I,90);
    end
    tmp = padarray(tmp,[1 1],0,'pre');
    
    I = padarray(I,[1 1],0,'pre');
    
    
    out = flattenMaskOverlay(bindVec(I),logical(tmp));
    
    
    imshow(out,[]);
    drawnow
    Y = [Y;tmp(:)];
    
    
    
    drawnow
end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Color images - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract images - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% extract bug eye
for e = 1:30
    fprintf(['start extracting:' num2str(e) '\n']);tic
    [fft{e}] = extractSingleBugEye_v2(setaria_FileList{e},15,[1 40]);
    fprintf(['done extracting:' num2str(e) ':' num2str(toc) '\n'])
end
%% stack data
imgTOT = 433*433;
RAD = size(fft{1}.f,3);
TH = size(fft{1}.f,4);
X = zeros([RAD TH 3 imgTOT*numel(fft)]);
str = 1;
for e = 1:numel(fft)
    
    fprintf(['start stacking:' num2str(e) ':' num2str(numel(fft)) '\n']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp = permute(fft{e}.f,[3 4 1 2]);
    tsz = size(tmp);
    tmp = reshape(tmp,[tsz(1:2) 433*433]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpA = fft{e}.A;
    tmpA = diff(cat(4,tmpA,tmpA(:,:,:,end)),1,4);
    tmpA = abs(tmpA)/(2*pi);
    tmpA = .5*(((tmpA.^2)+1).^.5 + (((1-tmpA).^2)+1).^.5);
    tmpA = permute(tmpA,[3 4 1 2]);
    tmpA = reshape(tmpA,[tsz(1:2) 433*433]);
   
    
    
    stp = str + size(tmpA,3) - 1;
    X(:,:,:,str:stp) = cat(3,reshape(tmp,[tsz(1:2) 1 prod(tsz(3:4))]),...
                             reshape(tmpA,[tsz(1:2) 1 prod(tsz(3:4))]),...
                             zeros([tsz(1:2) 1 prod(tsz(3:4))]));
    
    str =  stp + 1;
    fprintf(['end stacking:' num2str(e) ':' num2str(numel(fft)) '\n'])
end
%% convert to single
X = single(X);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract images - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Traing CNN - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define CNN
layers = [imageInputLayer([size(X,1) size(X,2) 3])
          convolution2dLayer([3,3],10)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          convolution2dLayer([3,3],10)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
      %{
layers = [imageInputLayer([size(X,1) size(X,2) 3])
          convolution2dLayer([7,2],10)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
      %}
NTR = round(1*numel(Y));
options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3);
%% train CNN
convnet_Setaria = trainNetwork(X(:,:,:,1:NTR),categorical(Y(1:NTR)),layers,options);
%save('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/setaria_CNN.mat','convnet_Setaria');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Traing CNN - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Run local and condor - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather data from CyVerse
% this is the RILs only path
dataPath = '/iplant/home/leakey_cyverse/setariaData/stomataTopoData%';  
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
setaria_FileList = {};
FileExt = {'nms'};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        setaria_FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%% search for file
for e = 1:numel(sorghum_FileList)
    if ~isempty(strfind(setaria_FileList{e},'601008 leaf2-2'))
        e
    end 
end
%% issue tickets over the FileList
[setaria_FileList] = issueBulkTicket(setaria_FileList);
%% issue ticket over the return folder
remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_setaria_ver0/'];
%remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_forCharles/'];
[remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numel(setaria_FileList),'write');
%% run test local
stomata_cnnVersion(setaria_FileList{100},convnet_Setaria,15,[1 40],'','');
%% launch on condor
func = cFlow('stomata_cnnVersion');
func.setMCRversion('v920');
func.setMemory('8000');
toRun = numel(setaria_FileList);
toRun = 300;
for e = 1:toRun
    func(setaria_FileList{e},convnet_Setaria,15,[1 40],'./output/',remoteOutputLocation);
    fprintf(['done rendering job:' num2str(e) ':' num2str(numel(setaria_FileList)) '\n'])
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
    
func.submitDag(auth,150,150);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Run local and condor - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setaria PIPELINE - end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% gather data from CyVerse
% this is the RILs only path
dataPath = '/iplant/home/leakey_cyverse/sorghumData%'; 
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
sorghum_FileList = {};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        sorghum_FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
dataPath = '/iplant/home/leakey_cyverse/maizeData%'; 
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
maize_FileList = {};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        maize_FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
dataPath = '/iplant/home/leakey_cyverse/setariaData%'; 
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
setariaData_FileList = {};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        setariaData_FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%%
CL{1} = sorghum_FileList;
CL{2} = maize_FileList;
CL{3} = setariaData_FileList;
mS = zeros(512,512,1,300);
cnt =1;
for e = 1:numel(CL)
    rndIDX = randperm(numel(CL{e}));
    
    for img = 1:100
        tmp = imread(CL{e}{rndIDX(img)});
        mS(:,:,1,cnt) = tmp;
        cnt = cnt + 1;
        img
    end
    e
end
YT = [ones(100,1);2*ones(100,1);3*ones(100,1)];

%%
rsM = [];
for e = 1:size(mS,4)
    rsM(:,:,:,e) = imresize(mS(:,:,:,e),.5);
end
%%
layers = [imageInputLayer([size(rsM,1) size(rsM,2) 1])
          convolution2dLayer([30,30],70)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          
          convolution2dLayer([10,10],30)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          
          convolution2dLayer([5,5],10)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          fullyConnectedLayer(3)
          softmaxLayer
          classificationLayer()];
 
      
NTR = round(numel(YT));
options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3);
%%
convnet_CL = trainNetwork(rsM,categorical(YT),layers,options);

%% my checks

%% scan for maize NMS files - RILS
FilePath1 = '/home/nate/Downloads/quickReturn_sorghum2016_ver1/';
iFileList1 = {};
FileExt = {'jpg'};
iFileList1 = gdig(FilePath1,iFileList1,FileExt,1);
FilePath2 = '/home/nate/Downloads/quickReturn_sorghum2016_ver2/';
iFileList2 = {};
FileExt = {'jpg'};
iFileList2 = gdig(FilePath2,iFileList2,FileExt,1);
%%
close all
clear nm1 nm2
for e = 1:numel(iFileList2)
    [p2,nm2{e}] = fileparts(iFileList2{e});
end

for e = 1:numel(iFileList2)
    [p1,nm1{e}] = fileparts(iFileList1{e});
end
[U,idx1,idx2] = intersect(nm1,nm2);
for u = 1:numel(U)
    I1 = imread(iFileList1{idx1(u)});
    I2 = imread(iFileList2{idx2(u)});
    imshow([I1,I2],[]);
    waitforbuttonpress
end
%% 
for e = 1:numel(iFileList)
    if isempty(strfind(iFileList{e},'_labeledImage'))
        rm(e) = true;
    else
        rm(e) = false;
    end
end
iFileList(rm) = [];
%%
h1 = figure;
h2 = figure;
for e = 1:numel(iFileList)
    tmp = imread(iFileList{e});
    cn1 = strrep(iFileList{e},'labeledImage.jpg','locations.csv');
    d = csvread(cn1);
    cnts(e,1) = size(d,1);
    figure(h2)
    [fp{e}] = impixel(tmp);
    [fn{e}] = impixel(tmp);
    cnts(e,2) = cnts(e,1) - size(fp{e},1) + size(fn{e},1);
    figure(h1);
    plot(cnts(:,2),cnts(:,1),'.')
    title(num2str(corr(cnts)))
end
%%
for e = 1:numel(fp)
    FF(e,:) = [size(fp{e},1) size(fn{e},1)];
end

            
