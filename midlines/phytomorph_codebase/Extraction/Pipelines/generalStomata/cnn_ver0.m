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
%load('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/generalStomata/RIL_sorghum_stomata_centers.mat','sorghumStomataCenter_row','sorghumStomataCenter_column');
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
    
    
    tmp2 = imdilate(tmp,strel('disk',11));
    tmp = imdilate(tmp,strel('disk',5,0));
    tmp = tmp2 + tmp;
    tmp = tmp2;
    
    for rot = 1:4
        tmp(1:border,:) = [];
        I(1:border,:) = [];
        tmp = imrotate(tmp,90);
        I = imrotate(I,90);
    end
    %tmp = padarray(tmp,[1 1],0,'pre');
    
    %I = padarray(I,[1 1],0,'pre');
    
    
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
    [fft{e}] = extractSingleBugEye_v2(sorghum_FileList{e},15,[3 20]);
    fprintf(['done extracting:' num2str(e) ':' num2str(toc) '\n'])
end
%% extract im2col style
I = imread(sorghum_FileList{1});
I = im2col(I,[41 41],'sliding');
X = zeros(41,41,1,size(I,2)*30);
Z = size(I,2);
str = 1;
for e = 1:30
    fprintf(['start extracting:' num2str(e) '\n']);tic
    I = imread(sorghum_FileList{e});
    I = im2col(I,[41 41],'sliding');
    I = reshape(I,[41 41 1 size(I,2)]);
    stp = str + Z - 1;
    X(:,:,:,str:stp) = I;
    str = stp + 1;
    fprintf(['done extracting:' num2str(e) ':' num2str(toc) '\n'])
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
          fullyConnectedLayer(3)
          softmaxLayer
          classificationLayer()];
%{
layers = [imageInputLayer([size(X,1) size(X,2) 3])
          convolution2dLayer([7,2],20)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
%}
      
NTR = round(1*numel(Y));
options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3);
%% diff rot invar
layers = [imageInputLayer([size(X,1) size(X,2) 1],'Normalization','None')
          convolution2dLayer([11,11],20)
          reluLayer
          maxPooling2dLayer([2 2],'Stride',2)
          convolution2dLayer([7,7],3)
          fullyConnectedLayer(2)
          softmaxLayer
          classificationLayer()];
imageAugmenter = imageDataAugmenter('RandRotation',[-90 90],'RandXScale',[.8 1.2],'RandYScale',[.8 1.2]);
%imageAugmenter = imageDataAugmenter('RandRotation',[-90 90]);
imageSize = [41 41 1];
NTR = 4;
TR = 25;
VAL = (30 - TR);
strTR = 1;
dSZ = 472;
stpTR = dSZ*dSZ*TR;
strVAL = stpTR + 1;
stpVAL = numel(Y);
VALData = {X(:,:,:,strVAL:NTR:stpVAL),categorical(Y(strVAL:NTR:stpVAL))};
datasource = augmentedImageSource(imageSize,X(:,:,:,strTR:NTR:end),categorical(Y(strTR:NTR:end)),'DataAugmentation',imageAugmenter);
%options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',3,'CheckpointPath','/mnt/snapper/nate/CP/','Plots','training-progress','ValidationFrequency',200,'ValidationData',VALData);
options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',8,'CheckpointPath','/mnt/snapper/nate/CP/','Plots','training-progress');
convnet_Sorghum6 = trainNetwork(datasource,net.Layers,options);
%convnet_Sorghum6 = trainNetwork(datasource,layers,options);
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
    if ~isempty(strfind(sorghum_FileList{e},'EF0366_1_1'))
        e
    end 
end
%% issue tickets over the FileList
[sorghum_FileList] = issueBulkTicket(sorghum_FileList);
%% issue ticket over the return folder
remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_sorghum2016_verFinal3/'];
%remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn_forCharles/'];
[remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numel(sorghum_FileList),'write');
%% run test local
figure
stomata_cnnVersion(sorghum_FileList{3521},convnet_Sorghum5,15,[3 20],'','');
%% run test local
figure
stomata_cnnVersion2(sorghum_FileList{3521},convnet_Sorghum6,15,[3 20],'','');
%% launch on condor
func = cFlow('stomata_cnnVersion');
func.setMCRversion('v920');
func.setMemory('12000');

toRun = numel(sorghum_FileList);
%toRun = 300;
for e = 1:toRun
    func(sorghum_FileList{e},convnet_Sorghum4,15,[3 40],'./output/',remoteOutputLocation);
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

            
