clear all

%% gather mat files from the run(s) for both cali and tera
% #1
FilePath = {};

FilePath{1} = '/mnt/tetra/nate/MN_MAT_LIST/All_030218-2018-03-02-20-00-22.9/';
FilePath{2} = '/mnt/tetra/nate/caliSample/';
for fl = 1:numel(FilePath)
    FileList{fl} = {};
    FileExt = {'mat'};
    FileList{fl} = gdig(FilePath{fl},FileList{fl},FileExt,1);
end
%% gather NEF files
% #2
imageFilePath = {};
imageFilePath{1} = '/mnt/tetra/nate/hirschSampleRAWWHOLE/';
imageFilePath{2} = '/mnt/tetra/nate/caliSampleRAWWHOLE/';
imageFileList = {};
for fl = 1:numel(imageFilePath)
    imageFileList{fl} = {};
    imageFileExt = {'nef'};
    imageFileList{fl} = gdig(imageFilePath{fl},imageFileList{fl},imageFileExt,1);
    imageFileList{fl} = imageFileList{fl}(randperm(numel(imageFileList{fl})));
end
FileList = imageFileList;
%% make whole list
% #3
cnt = 1;
for fl = 1:numel(imageFileList)
    for e = 1:numel(imageFileList{fl})
        wholeList{cnt} = imageFileList{fl}{e};
        cnt = cnt + 1;
    end
end
%% 
metaTable = myDS();
metaTable.importImagesFromList(wholeList,'maizeSeedling_rawImage');
%% get list from irods
dataPath = '/iplant/home/hirsc213/maizeData/seedlingData/';
FileList = myDS.igdig_irods(dataPath);
%% import from irods
metaTable = myDS('/mnt/tetra/nate/seedlingImageParts/');
tic
metaTable.importImagesFromList(FileList(1:1000),'maizeSeedling_rawImage');
toc
%% TEST for reporting string
metaTable.toString(true);
%% TEST insert queue
qKey = metaTable.initQueue('string');
metaTable.isQueueEmpty(qKey);
metaTable.pushToDataQueue(qKey,'quedata1','string');
metaTable.pushToDataQueue(qKey,'quedata2','string');
metaTable.pushToDataQueue(qKey,'quedata3','string');
%% create TEST context
[contextKey] = metaTable.initExecuteContext(1);
%%
fracDpi = 1;
 % set to default value of 10^6
defaultAreaPix = 10^6;
defaultAreaPix = round(defaultAreaPix*fracDpi);
% set to default value of 300
rho = 300;
rho = round(rho*fracDpi);
% set to default value of 70
colRange1 = 70;
% set to default value of 166
colRange2 = 166;
% set to default value of 50
fill = 50;
toSave = false;
toDisplay = false;
func = @(X)singleCobImage(X,numberOfObjects,oPath,remotePath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay);

metaTable.initAlgoContext(func);

%%
metaTable.initExecuteSequenceByType('maizeSeedling_rawImage','maizePipeline');
%%
sub = metaTable.partition(3,1);
%%
IDX = findImageFramePair(metaTable,'maizeSeedling_rawImage');
%%
baseLocation = '/mnt/tetra/nate/seedlingImageParts/';
tT = {};
for e = 1:5%0%:size(IDX,1)
    try
        tT{e} = generateSceneData_ver2(sub,obj,nonBackgroundGMM,125,20,400,baseLocation);
        e
    catch
    end
end

%% loop load over the mat files over the groups
% #4
N = 200;
%totD = 100;
n = 1;
parfor fl = 1:numel(FileList)
    totD = numel(FileList{fl});
    totD = 800;
    %totD = 1;
    KS{fl} = zeros(3*N*totD,3);
    MS{fl} = zeros(3*N*totD,1);
    str = 1;
    stp = str + N - 1;
    
    for e = 1:totD
        try
            tmp = load(FileList{fl}{e});
            % for each plant
            for i = 1:numel(tmp.returnI)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % look at the size of the crop box
                sz = size(tmp.returnI{i});
                % reshape the data
                tmp1 = reshape(tmp.returnI{i},[prod(sz(1:2)) sz(3)]);
                tmp2 = reshape(tmp.MASK{i},[prod(sz(1:2)) 1]);
                UQ = unique(tmp2);
                % for each labled type
                
                
                for u = 1:numel(UQ)
                    % find the pixels labled as such
                    fidx = find(tmp2 == UQ(u));
                    % permute the labels
                    fidx = fidx(randperm(numel(fidx)));
                    KS{fl}(str:stp,:) = tmp1(fidx(1:N),:);
                    MS{fl}(str:stp,:) = tmp2(fidx(1:N),:);
                    
                    str = stp + 1;
                    stp = str + N - 1;
                end
                
               

            end
        catch ME
            ME
            %toRM{n} = str:stp;
            %n = n + 1;
        end
        e
    end
    
end
%% build up secondary distribution
% #4.5
colorSample_oPath = '/mnt/tetra/nate/phunSeedlings/colorSample/';
FilePath = colorSample_oPath;
cFileList = {};
FileExt = {'mat'};
cFileList = gdig(FilePath,cFileList,FileExt,1);
CS2 = [];
for e = 1:numel(cFileList)
    load(cFileList{e},'colorSample');
    rn = randperm(size(colorSample,1));
    CS2 = [CS2;colorSample(1:rn(3000),:)];
    e
end
%%
obj2nd = gmdistribution.fit(CS2(1:10:end,:),6);
%% view the data
% #5
close all
CL = {'b.' 'c.' 'r.' 'm.'};
cnt = 1;
masterColorStack = [];
masterColorLabel = [];
for e = 1:numel(KS)
    UQ = unique(MS{e});
    masterColorStack = [masterColorStack;KS{e}];
    masterColorLabel = [masterColorLabel;MS{e}];
    for u = 1:numel(UQ)
        fidx = find(MS{e}==UQ(u));
        plot3(KS{e}(fidx(1:100:end),1),KS{e}(fidx(1:100:end),2),KS{e}(fidx(1:100:end),3),CL{cnt});
        cnt = cnt +1;
        hold on
    end
end
%% cluster the data back ground er no
% #6
UQ = unique(masterColorLabel);
mu2 = [];
sigma2 = [];
GRPS = [2 1];
colorGroups = {};
for u = 1:numel(UQ)
    fidx = find(masterColorLabel == UQ(u));
    N(u) = numel(fidx);
    rm = all(masterColorStack(fidx,:) == 0);
    fidx(rm) = [];
    mu(u,:) = mean(masterColorStack(fidx,:));
    sigma(:,:,u) = cov(masterColorStack(fidx,:));
    
    GMModel{u} = fitgmdist(masterColorStack(fidx,:),GRPS(u));
    
    mu2 = [mu2;GMModel{u}.mu];
    sigma2 = cat(3,sigma2,GMModel{u}.Sigma);
    
    colorGroups{u} = masterColorStack(fidx,:);
end
p = N / sum(N);
obj = gmdistribution(mu,sigma,p);
%obj = gmdistribution(mu2,sigma2);
%% find background - #7
backGround_oPath = '/mnt/tetra/nate/phunSeedlings2/backgrounds/';
qrcode_oPath = '/mnt/tetra/nate/phunSeedlings2/qrcodes/';
colorSample_oPath = '/mnt/tetra/nate/phunSeedlings2/colorSample/';
extraData_oPath = '/mnt/tetra/nate/phunSeedlings2/extraTop/';
conetainerData_oPath = '/mnt/tetra/nate/phunSeedlings2/conetainer_EX/';
mkdir(backGround_oPath);
mkdir(qrcode_oPath);
mkdir(colorSample_oPath);
mkdir(extraData_oPath);
mkdir(conetainerData_oPath);
cnt = 1;
close all
for fl = 1:2
    parfor e = 1:800
        try
            nm = (fl-1)*400 + e;
            stripComponents1(imageFileList{fl}{e},obj,obj2nd,num2str(nm),backGround_oPath,qrcode_oPath,colorSample_oPath,extraData_oPath,conetainerData_oPath,false);
            cnt = cnt + 1;
        catch
        end
    end
end
%% GENERATE setup

imgSZ = [4020 6036];
oSIM = '/mnt/tetra/nate/phunSeedlings2/sims2/';
mkdir(oSIM);
backGround_oPath = '/mnt/tetra/nate/phunSeedlings2/backgrounds/';
qrcode_oPath = '/mnt/tetra/nate/phunSeedlings2/qrcodes/';
colorSample_oPath = '/mnt/tetra/nate/phunSeedlings2/colorSample/';
extraData_oPath = '/mnt/tetra/nate/phunSeedlings2/extraTop/';
conetainerData_oPath = '/mnt/tetra/nate/phunSeedlings2/conetainer/';
conetainerData_oPath = '/mnt/tetra/nate/phunSeedlings2/SIMconetainer/';


FilePath = backGround_oPath;
bkFileList = {};
FileExt = {'mat'};
bkFileList = sdig(FilePath,bkFileList,FileExt,1);
close all

FilePath = qrcode_oPath;
qrFileList = {};
FileExt = {'mat'};
qrFileList = gdig(FilePath,qrFileList,FileExt,1);


FilePath = conetainerData_oPath;
cnFileList = {};
FileExt = {'mat'};
cnFileList = gdig(FilePath,cnFileList,FileExt,1);
idx1 = contains(cnFileList,'_1');
idx2 = contains(cnFileList,'_2');
idx3 = contains(cnFileList,'_3');
coneFileList{1} = cnFileList(idx1);
coneFileList{2} = cnFileList(idx2);
coneFileList{3} = cnFileList(idx3);

close all
%% run GENERATE
parfor R = 1:10000
    simScene(bkFileList,qrFileList,coneFileList,imgSZ,oSIM,R,false);
    R
end
%% load data for cone type network
TEST = true;
if ~TEST
    coneI = [];
end
imgSZ = [4020 6036];
for e = 1:5000
    rnd = randi(numel(cnFileList),1);
    q1 = randi(numel(qrFileList),1);
    d = load(cnFileList{rnd});
    d.coneMask = imclose(d.coneMask,strel('square',121));
    d.coneMask = imfill(d.coneMask,'holes');
    
    
        QR = load(qrFileList{q1},'QR','boundingBox','MQR');
        QRLOC(1) = randi(imgSZ(2) - 500);
        QRLOC(2) = randi(round(randi(imgSZ(1)/2)));
        QR.boundingBox(1:2) = QRLOC;

        n1 = randi(numel(bkFileList{1}));
        n2 = randi(numel(bkFileList{2}));
        n3 = randi(numel(bkFileList{3}));
        n4 = randi(numel(bkFileList{4}));
        
        load(bkFileList{1}{n1},'hB1')
        load(bkFileList{2}{n2},'hB2')
        load(bkFileList{3}{n3},'vB1')
        load(bkFileList{4}{n4},'vB2')

        [simBK] = generateBackground(vB1,hB1,vB2,hB2,imgSZ);
    
        
        
        TOP_blend = zeros(size(simBK,1),size(simBK,2));
        TOP_color = zeros(size(simBK));
        
    
        d.toDropMask = d.coneMask;
        [TOP_color,TOP_blend] = dropImage(TOP_color,d.boundingBox,d.coneImage,d.toDropMask,TOP_blend);

        
        TOP_blend_blur = imfilter(TOP_blend,fspecial('disk',11),'replicate');
        TOP_blend = TOP_blend.*TOP_blend_blur;
        BOT_blend = 1 - TOP_blend;


        TOP_color = bsxfun(@times,TOP_color,TOP_blend);
        BOT_color = bsxfun(@times,simBK,BOT_blend);
        TOT_color = BOT_color + TOP_color;
        
        newR = regionprops(logical(TOP_blend),'BoundingBox');
        tmp = imcrop(TOT_color,newR(1).BoundingBox);
    W = tmp;
    tmp = imresize(tmp,[1000 1500]);
    tmp = imresize(tmp,.05);
    if TEST 
        cT = coneTypeNetwork.classify(tmp);
        imshow(W,[]);
        title(num2str(double(cT)))
        drawnow
    else
         coneI(:,:,:,e) = tmp;
        isCone(e) = contains(cnFileList{rnd},'without');
         imshow(tmp,[]);
    drawnow
    end
   
   
    e
end
%% learn cone type
imgSZ = size(coneI);
imgSZ(4) = [];
inputLayer = imageInputLayer([imgSZ]);
middleLayer1 = [...
    convolution2dLayer([10 10],12)
    reluLayer
    maxPooling2dLayer([5 5],'Stride',2)];
middleLayer2 = [...
    convolution2dLayer([5 5],4)
    reluLayer
    maxPooling2dLayer([2 2],'Stride',5)];
middleLayer3 = [...
    convolution2dLayer([5 1],3)
    reluLayer
    maxPooling2dLayer([2 1],'Stride',2)];
finalLayers = [...
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
layers = [...
    inputLayer
    middleLayer1
    middleLayer2
    finalLayers];
options = trainingOptions('sgdm',...
    'MaxEpochs',1000, ...
    'Verbose',true,...
    'InitialLearnRate',.05,...
    'Plots','training-progress',...
    'MiniBatchSize',128*1,...
    'ExecutionEnvironment','parallel');
coneTypeNetwork = trainNetwork(coneI,categorical(~logical(isCone')), layers, options);
%% build up secondary distribution
colorSample_oPath = '/mnt/tetra/nate/phunSeedlings/colorSample/';
FilePath = colorSample_oPath;
cFileList = {};
FileExt = {'mat'};
cFileList = gdig(FilePath,cFileList,FileExt,1);
CS2 = [];
for e = 1:numel(cFileList)
    load(cFileList{e},'colorSample');
    rn = randperm(size(colorSample,1));
    CS2 = [CS2;colorSample(1:rn(3000),:)];
    e
end

%% LOAD QR SCANDLINE MODERN
simcode_oPath = '/mnt/tetra/nate/phunSeedlings/sims2/';
%simmini_oPath = '/mnt/tetra/nate/phunSeedlings/sims_mini/';
simFilePath = simcode_oPath;
simFileList = {};
FileExt = {'tif'};
simFileList = gdig(simFilePath,simFileList,FileExt,1);
Y = [];
X = [];

xBOXH = [];
yBOXH = [];

xBOXV = [];
yBOXV = [];

imgSZ = [4020 6036];
trainTable = table;
ptr = 1;
scanlineSZ = 5;

toLoad = 800;

I = imread(simFileList{ptr});
I = double(imresize(I,.03))/255;
matFile = strrep(simFileList{ptr},'.tif','.mat');
d = load(matFile);
tmpMask = d.wholeQRMask(1:imgSZ(1),1:imgSZ(2));
tmpMask = imresize(tmpMask,.03);

coneMask = d.TOP_blend(1:imgSZ(1),1:imgSZ(2));
coneMask = imclose(logical(coneMask),strel('disk',21,0));
bot = coneMask(end,:);
coneMask(end,:) = 1;
coneMask = imfill(coneMask,'holes');
coneMask(end,:) = bot;
coneMask = imresize(coneMask,.03,'nearest');
coneMask = edge(coneMask,'Sobel',[],'vertical');
coneMask = repmat(coneMask(end,:),[size(coneMask,1) 1]);
coneMask = double(coneMask);
CscanlineSZ = 13;
[tmp_CLINE_Y,tmp_CLINE_X,~,~] = generateScanlineData(I((end-1)/2:end,:,:),coneMask((end-1)/2:end,:),CscanlineSZ,size(I((end-1)/2:end,:,:)));
SKIPCY = size(tmp_CLINE_Y,1);
xBOXC = zeros([size(tmp_CLINE_X,1) size(tmp_CLINE_X,2) size(tmp_CLINE_X,3) SKIPCY*toLoad]);
yBOXC = zeros(SKIPCY*toLoad,size(tmp_CLINE_Y,2));

[tmp_HLINE_Y,tmp_HLINE_X,tmp_VLINE_Y,tmp_VLINE_X] = generateScanlineData(I,tmpMask,scanlineSZ,size(I));
SKIPHY = size(tmp_HLINE_Y,1);
SKIPVY = size(tmp_VLINE_Y,1);
yBOXH = zeros(SKIPHY*toLoad,size(tmp_HLINE_Y,2));
yBOXV = zeros(SKIPVY*toLoad,size(tmp_VLINE_Y,2));
xBOXH = zeros([size(tmp_HLINE_X,1) size(tmp_HLINE_X,2) size(tmp_HLINE_X,3) SKIPHY*toLoad]);
xBOXV = zeros([size(tmp_VLINE_X,1) size(tmp_VLINE_X,2) size(tmp_VLINE_X,3) SKIPVY*toLoad]);
strH = 1;
strV = 1;
strC = 1;
numL = 1;
%for e = 1:1000%00%numel(simFileList)
ptr = 1;
while numL <= toLoad
    try
        
        I = imread(simFileList{ptr});
        I = double(imresize(I,.03))/255;
        matFile = strrep(simFileList{ptr},'.tif','.mat');
        d = load(matFile);

        if d.toQR
            
            stpH = strH + SKIPHY - 1;
            stpV = strV + SKIPVY - 1;
            stpC = strC + SKIPCY - 1;

            tmpMask = d.wholeQRMask(1:imgSZ(1),1:imgSZ(2));
            tmpMask = imresize(tmpMask,.03);
            
            coneMask = d.TOP_blend(1:imgSZ(1),1:imgSZ(2));
            coneMask = imclose(logical(coneMask),strel('disk',21,0));
            bot = coneMask(end,:);
            coneMask(end,:) = 1;
            coneMask = imfill(coneMask,'holes');
            coneMask(end,:) = bot;
            coneMask = imresize(coneMask,.03,'nearest');
            coneMask = edge(coneMask,'Sobel',[],'vertical');
            coneMask = repmat(coneMask(end,:),[size(coneMask,1) 1]);
            coneMask = double(coneMask);
            [tmp_CLINE_Y,tmp_CLINE_X,~,~] = generateScanlineData(I((end-1)/2:end,:,:),coneMask((end-1)/2:end,:),CscanlineSZ,size(I((end-1)/2:end,:,:)));
            yBOXC(strC:stpC,:) = tmp_CLINE_Y;
            xBOXC(:,:,:,strC:stpC) = tmp_CLINE_X;
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [tmp_HLINE_Y,tmp_HLINE_X,tmp_VLINE_Y,tmp_VLINE_X] = generateScanlineData(I,tmpMask,scanlineSZ,size(I));
            yBOXH(strH:stpH,:) = tmp_HLINE_Y;
            xBOXH(:,:,:,strH:stpH) = tmp_HLINE_X;
            yBOXV(strV:stpV,:) = tmp_VLINE_Y;
            xBOXV(:,:,:,strV:stpV) = tmp_VLINE_X;
            %{
            yBOXH = [yBOXH;tmp_HLINE_Y];
            xBOXH = cat(4,xBOXH,tmp_HLINE_X);
            yBOXV = [yBOXV;tmp_VLINE_Y];
            xBOXV = cat(4,xBOXV,tmp_VLINE_X);
            %}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            strH = stpH + 1;
            strV = stpV + 1;
            strC = stpC + 1;
            numL = numL + 1;
            numL
        end
        Y(ptr,1) = d.toQR;
        Y(ptr,2) = d.numCONE;
        Y(ptr,3:4) = d.countType;
        X(:,:,:,ptr) = I;
        imshow(I,[]);
        title([num2str(ptr) '-->' num2str(numL)]);
        drawnow
        ptr = ptr + 1;
    catch ME
        ptr = ptr + 1;
        ME
    end
end
%% learn box part 1 for QR code
imgSZ = size(xBOXH);
imgSZ(4) = [];
inputLayer = imageInputLayer([imgSZ]);
middleLayer1 = [...
    convolution2dLayer([10 5],12)
    reluLayer
    maxPooling2dLayer([5 1],'Stride',2)];
middleLayer2 = [...
    convolution2dLayer([5 1],3)
    reluLayer
    maxPooling2dLayer([2 1],'Stride',2)];
middleLayer3 = [...
    convolution2dLayer([5 1],3)
    reluLayer
    maxPooling2dLayer([2 1],'Stride',2)];
finalLayers = [...
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
layers = [...
    inputLayer
    middleLayer1
    finalLayers];
options = trainingOptions('sgdm',...
    'MaxEpochs',100, ...
    'Verbose',true,...
    'InitialLearnRate',.05,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
QRboxH = trainNetwork(xBOXH,categorical(yBOXH), layers, options);
%% TRAIN box V net
imgSZ = size(xBOXV);
imgSZ(4) = [];
inputLayer = imageInputLayer([imgSZ]);
layers = [...
    inputLayer
    middleLayer1
    finalLayers];
QRboxV = trainNetwork(xBOXV,categorical(yBOXV), layers, options);
%% learn box part 1 for QR code
imgSZ = size(xBOXC);
imgSZ(4) = [];
inputLayer = imageInputLayer([imgSZ]);
middleLayer1 = [...
    convolution2dLayer([15 13],12)
    reluLayer
    maxPooling2dLayer([5 1],'Stride',2)];
middleLayer2 = [...
    convolution2dLayer([5 1],3)
    reluLayer
    maxPooling2dLayer([2 1],'Stride',2)];
middleLayer3 = [...
    convolution2dLayer([5 1],3)
    reluLayer
    maxPooling2dLayer([2 1],'Stride',2)];
finalLayers = [...
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
layers = [...
    inputLayer
    middleLayer1
    middleLayer2
    finalLayers];

options = trainingOptions('sgdm',...
    'MaxEpochs',1000, ...
    'Verbose',true,...
    'InitialLearnRate',.05,...
    'Plots','training-progress',...
    'MiniBatchSize',128*4,...
    'ExecutionEnvironment','parallel');
z1 = sum(yBOXC);
f1 = find(yBOXC==1);
f0 = find(yBOXC==0);
f0 = f0(randperm(numel(f0)));
N = 2;
fM = [f1;f0(1:N*z1)];
QRboxC = trainNetwork(xBOXC(:,:,:,fM),categorical(yBOXC(fM,:)), layers, options);
%% QR detection
imgSZ = size(X);
imgSZ(4) = [];
inputLayer = imageInputLayer([imgSZ]);
middleLayer1 = [...
    convolution2dLayer([20 20],12)
    reluLayer
    crossChannelNormalizationLayer(4)
    maxPooling2dLayer([5 5],'Stride',10)];
middleLayer2 = [...
    convolution2dLayer(5,3)
    reluLayer
    crossChannelNormalizationLayer(4)
    maxPooling2dLayer(2,'Stride',2)];
finalLayers = [...
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
layers = [...
    inputLayer
    middleLayer1
    middleLayer2
    finalLayers];
options = trainingOptions('sgdm',...
    'MaxEpochs',100, ...
    'Verbose',true,...
    'InitialLearnRate',.05,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
%% TRAIN QR network - BOOLEAN
toSKIP = 5;
trainedNet = trainNetwork(X,categorical(Y(:,1)),layers,options);
%% init count network
middleLayer1 = [...
    convolution2dLayer([20 20],3)
    reluLayer
    crossChannelNormalizationLayer(4)
    maxPooling2dLayer([5 5],'Stride',10)];
finalLayersR = [...
    fullyConnectedLayer(3)
    regressionLayer];
Rlayers = [...
    inputLayer
    middleLayer1
    finalLayersR];

options = trainingOptions('sgdm',...
    'MaxEpochs',100, ...
    'Verbose',true,...
    'InitialLearnRate',.01,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
%% TRAIN count network 
toSKIP = 1;
trainedNetR2 = trainNetwork(X(:,:,:,1:toSKIP:end),Y(1:toSKIP:end,2:end),Rlayers,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test networks
close all
simcode_oPath = '/mnt/tetra/nate/phunSeedlings/sims2/';
simFilePath = simcode_oPath;
testFileList = {};
FileExt = {'tif'};

testFileList = gdig(simFilePath,testFileList,FileExt,1);

%%
close all
FilePath = '/mnt/tetra/nate/seedlingDATApile/maizeData/seedlingData/';
FileList = {};
FileExt = {'nef'};
testFileList = gdig(FilePath,FileList,FileExt,1);

Ypre = [];
Ypost = [];

    toConnect{1} = [[3 5];[5 7];[7 8];[8 3]];
    toConnect{2} = [[5 1];[1 6];[6 7];[7 5]];
    toConnect{3} = [[8 6];[6 2];[2 4];[4 8]];
%testFileList = wholeFileList;

qrPilePath = '/mnt/tetra/nate/maizeSeedlingData/qrPile/';
thumbPilePath = '/mnt/tetra/nate/maizeSeedlingData/thumb_image_pile/';
tableFile = '/mnt/tetra/nate/maizeSeedlingData/tableFile.mat';
trainTable = table;

mkdir(qrPilePath);
mkdir(thumbPilePath);
TOT =  1;





for e = 1:numel(testFileList)
    try
        I = double(imread(testFileList{e}))/255;
        I = imresize(I,imgSZ);
        Io = I;
        I = imresize(I,.03);
        matFile = strrep(testFileList{e},'.tif','.mat');
        %d = load(matFile);
        tic
        coneCount = round(trainedNetR2.predict(I));
        toc
        tic
        QRis = double((trainedNet.classify(I)))-1;
        toc
        
        [~,tmp_CLINE_X,~,~] = generateScanlineData(I((end-1)/2:end,:,:),[],CscanlineSZ,size(I((end-1)/2:end,:,:)));
           
        
        
        IB = imresize(double(Io),.03);
        [tmp_HLINE_Y,tmp_HLINE_X,tmp_VLINE_Y,tmp_VLINE_X] = generateScanlineData(IB,[],scanlineSZ,size(IB));
        byX = double(QRboxH.classify(tmp_HLINE_X))-1;
        byY = double(QRboxV.classify(tmp_VLINE_X))-1;
        BY = byY(:,1)*byX(:,1)';
        BY = padarray(BY,[2 2],0,'both');
        BY = imresize(BY,[size(Io,1) size(Io,2)]);
        R = regionprops(BY > .5^2,'BoundingBox','Area');
        R([R.Area] < 700000) = [];
        
        
        
        SIDEWALLS = QRboxC.classify(tmp_CLINE_X);
        SR = regionprops(logical(double(SIDEWALLS)-1));
        swl = [];
        for sw = 1:numel(SR)
            swl(sw,:) = SR(sw).Centroid;
            swl(sw,2) = swl(sw,2) + (CscanlineSZ-1)/2-1;
        end
        swl(:,2) = swl(:,2)*.03^-1;
            
        
        ppp = [];
        if ~isempty(R)
            
             %[dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(oI,20);
            
            
            
            QR_crop = imcrop(Io,R(1).BoundingBox);
           
           
            
            
            oSZ = size(QR_crop);
            QR_cropR = imresize(QR_crop,qrNORSZ);
            toScale = size(QR_cropR).*oSZ.^-1;
            oSZ2 = size(QR_cropR);
            QR_cropR = imresize(QR_cropR,qrSZ);
            toScale2 = size(QR_cropR).*oSZ2.^-1;
            
            [markX,tmp_squareY] = generateSquarescanData(QR_cropR,[],patchSZ,[]);
            [MARKS,MARKS_V] = trainedNet_QR_mark2.classify(markX,'MiniBatchSize',256*512);
            V = reshape(double(MARKS(:,1)),[size(QR_cropR,1) size(QR_cropR,2)]-(patchSZ-1));
            V_V = reshape(double(MARKS_V),[[size(QR_cropR,1) size(QR_cropR,2)]-(patchSZ-1) size(MARKS_V,2)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % gather the points for the red border
            ppp = [];
            red_corner_save = [];
            for pp = 1:8
                tmpV = V == pp;
                tmpV = imopen(tmpV,strel('disk',2,0));
                tmpR = regionprops(tmpV,'Centroid','Area');
                [~,sidx] = sort([tmpR.Area],'descend');
                tmpR = tmpR(sidx);
                ppp(pp,:) = tmpR(1).Centroid;
            end
            ppp = bsxfun(@plus,ppp,(patchSZ-1)/2);
            %{
            figure;
            imshow(QR_cropR,[]);
            hold on
            plot(ppp(:,1),ppp(:,2),'r*')
            %}
            ppp = bsxfun(@times,ppp,fliplr(toScale2(1:2)).^-1);
            ppp = bsxfun(@times,ppp,fliplr(toScale(1:2)).^-1);
            
            
            
            ppp = bsxfun(@plus,ppp,R(1).BoundingBox(1:2));
            my_corner_order = [3 8 4 5 7 1 6 2];
            red_corner_save = ppp(my_corner_order,:);
            %trainTable(TOT,'metaDataBox1') = {{R(1).BoundingBox}};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % gather the points for the red border
            pppB = [];
            blue_corner_save = [];
            cnt = 1;
            for pp = 9:12
                tmpV = V == pp;
                tmpV = imopen(tmpV,strel('disk',2,0));
                tmpR = regionprops(tmpV,V_V(:,:,pp),'Centroid','Area','MeanIntensity');
                [~,sidx] = sort([tmpR.MeanIntensity],'descend');
                tmpR = tmpR(sidx);
                for c = 1:24
                    pppB(c,:,cnt) = tmpR(c).Centroid;
                end
               
                
                blue_corner_order = [4 1 2 3];
                
                pppB(:,:,cnt) = bsxfun(@plus,pppB(:,:,cnt),(patchSZ-1)/2);
                pppB(:,:,cnt) = bsxfun(@times,pppB(:,:,cnt),fliplr(toScale2(1:2)).^-1);
                pppB(:,:,cnt) = bsxfun(@times,pppB(:,:,cnt),fliplr(toScale(1:2)).^-1);
                pppB(:,:,cnt) = bsxfun(@plus,pppB(:,:,cnt),R(1).BoundingBox(1:2));
                
                
                
                cnt = cnt + 1;
            end
            
            blue_corner_save = pppB(:,:,blue_corner_order);
            
            
            
            % HIGH RES BOXES
            fN = [thumbPilePath num2str(e) '.tif'];
            imwrite(Io,[thumbPilePath num2str(e) '.tif']);
            imwrite(Io,[thumbPilePath num2str(e) '.tif']);
            trainTable(TOT,'wholeImageFileName') = {fN};
           
            
            
            SNAP = 40;
            SNIP = 2.6;
            JERK = round(SNAP*SNIP);
            tightBB = [ppp(3,:)-JERK ppp(1,1)-ppp(3,1)+JERK*2 ppp(4,2)-ppp(3,2)+JERK*2];
            staticBB = [ppp(3,:)-JERK ppp(5,1)-ppp(3,1)+JERK*2 ppp(8,2)-ppp(3,2)+JERK*2];
            dynamicBB = [ppp(8,:)-JERK ppp(1,1)-ppp(8,1)+JERK*2 ppp(4,2)-ppp(8,2)+JERK*2];
            detectedImg = insertShape(Io, 'Rectangle', tightBB,'Color',{'m'},'LineWidth',8);
            detectedImg = insertShape(detectedImg, 'Rectangle', staticBB,'Color',{'g'},'LineWidth',8);
            detectedImg = insertShape(detectedImg, 'Rectangle', dynamicBB,'Color',{'g'},'LineWidth',8);
            
            
            
            for blah = 1:size(blue_corner_save,3)
                blue_corner_save(:,:,blah) = bsxfun(@minus,blue_corner_save(:,:,blah),tightBB(1:2));
            end
            
            red_corner_save = bsxfun(@minus,red_corner_save,tightBB(1:2));
            trainTable(TOT,'metaDataBox') = {{tightBB}};
            trainTable(TOT,'staticDataBox') = {{staticBB}};
            trainTable(TOT,'dynamicDataBox') = {{dynamicBB}};
            trainTable(TOT,'dynamicPoints') = {{blue_corner_save}};
            trainTable(TOT,'redPoints') = {{red_corner_save}};
            
            
            
            
            QR_crop_save = imcrop(Io,tightBB);
           
            
            imwrite(QR_crop_save,[qrPilePath num2str(e) '.tif']);
            trainTable(TOT,'qrImageFileName') = {[qrPilePath num2str(e) '.tif']};
            save(tableFile,'trainTable');
            
            
            
            
            TOT = TOT + 1;
            
            
            
            
        end
        
        
        
        
        
        Ypre(e,2:4) = coneCount;
        Ypre(e,1) = QRis;
        %Ypost(e,1) = double(d.toQR);
        %Ypost(e,2) = double(d.numCONE);
        %Ypost(e,3:4) = double(d.countType);
        %corr(Ypost,Ypre)
        %if Ypost(e,2) ~= Ypre(e,2)
            imshow(detectedImg,[]);
            title(['isQR:' num2str(QRis) ':' num2str(Ypre(e,2:4)) '-->' num2str(Ypre(e,2:4))]);
            %title(['isQR:' num2str(QRis) ':' num2str(Ypre(e,2:4)) '-->' num2str(Ypost(e,2:4))]);
            if ~isempty(R)
                rectangle('Position',R(1).BoundingBox);
            end
            if ~isempty(ppp)
                hold on
                plot(ppp(:,1),ppp(:,2),'r*')
                for bb = 1:numel(toConnect)
                    for l = 1:size(toConnect{bb},1)
                        plot(ppp(toConnect{bb}(l,:),1),ppp(toConnect{bb}(l,:),2),'k')
                    end
                end
                
                hold on
                plot(pppB(:,1,1),pppB(:,2,1),'b*')
                plot(pppB(:,1,2),pppB(:,2,2),'b*')
                plot(pppB(:,1,3),pppB(:,2,3),'b*')
                plot(pppB(:,1,4),pppB(:,2,4),'b*')
                
                
                
                hold off
            end
            
            hold on
            for sw = 1:size(swl,1)
                plot(swl(sw,2)*ones(size(Io,1),1),1:size(Io,1),'r');
            end
            hold off
            drawnow
            
        %end
    catch ME
        ME
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
obj2nd = gmdistribution.fit(CS2(1:10:end,:),6);
%% model QR code
qrcode_oPath = '/mnt/tetra/nate/phunSeedlings2/qrcodes/';
FilePath = qrcode_oPath;
QRFileList = {};
FileExt = {'mat'};
QRFileList = gdig(FilePath,QRFileList,FileExt,1);
%% mark QR sheet
close all
clear toConnect
sel = [2 1];
sel = [5 3 2];
sel = [1 3 5];
h1 = figure;
h2 = figure;
boxCNT = 1;
blueBOX_stack = [];
blueBOX_size = [];
LOC = [];
frameY = [];
frameX = [];
ref = load(QRFileList{1},'QR');
ref = ref.QR;
for e = 1:numel(QRFileList)
    try
        clear b PO skelP
        load(QRFileList{e},'QR')
        sz = size(QR);
        I = QR;
        
        %QR = imhistmatch(QR,ref);
        GI = rgb2gray(I);
        Cpoints = detectHarrisFeatures(GI);

        QR = reshape(QR,[prod(sz(1:2)) sz(3)]);
        [k,~,KK] = QRlabel2_2.cluster(QR);
        k = reshape(k,sz(1:2));
        KK = reshape(KK,[sz(1:2) size(KK,2)]);
        rgb = label2rgb(k);
        frameMSK = k == sel(1);
        frameMSK = frameMSK.*(KK(:,:,sel(1)) > .999);
        paperMSK = k ~= sel(2);
        blueSquareMSK = KK(:,:,sel(3)) > .98;
        blueSquareMSK = bwareaopen(blueSquareMSK,100);

        %blueSquareMSK = bwlarge(blueSquareMSK,24);
        frameMSK = bwareaopen(frameMSK,300);
        frameMSK = imclose(frameMSK,strel('line',20,0));
        frameMSK = imclose(frameMSK,strel('line',20,90));
        %frameMSK = imerode(frameMSK,strel('square',5));
        frameMSK = bwlarge(frameMSK);


        frameRegions = imfill(frameMSK,'holes') - frameMSK;
        frameRegions1 = logical(bwlarge(frameRegions,1));
        frameRegions2 = logical(bwlarge(frameRegions,2) - frameRegions1);
        frameRegions3 = logical(bwlarge(frameRegions,3) - frameRegions2 - frameRegions1);



        blueSquareMSK = logical(blueSquareMSK.*frameRegions1);

        blueSquareMSK = imfill(blueSquareMSK,'holes');
        blueSquareMSK = bwlarge(blueSquareMSK,24);
        blueSQUARE_REG = regionprops(blueSquareMSK,'PixelIdxList','BoundingBox','Centroid');

        BOXCEN = [];
        for box = 1:numel(blueSQUARE_REG)
            BOXCEN(box,:) = blueSQUARE_REG(box).Centroid;
            img = imcrop(I,blueSQUARE_REG(box).BoundingBox);
            blueBOX_size(boxCNT,:) = [size(img,1) size(img,2)];
            img = imresize(img,[65 80]);
            blueBOX_stack(:,:,:,boxCNT) = img;
            %figure(h2);
            %imshow(img,[]);
            %boxCNT = boxCNT + 1;
            %drawnow
        end

        nBOXSTACK = [];
        for col = 1:6
            [~,bidx] = sort(BOXCEN(:,1));
            [~,bidx2] = sort(BOXCEN(bidx(1:4),2));
            nBOXSTACK = [nBOXSTACK;BOXCEN(bidx(bidx2),:)];
            BOXCEN(bidx(bidx2),:) = [];
        end
        BOXCEN = nBOXSTACK;


        paperMSK = bwlarge(paperMSK);
        paperMSK = imfill(paperMSK,'holes');
        paperR = regionprops(paperMSK);


        innerMSK = imfill(frameMSK,'holes');


        LIGHTING = k == 4;

        HSV = rgb2hsv(I);
        GRAD = imfilter(HSV(:,:,3),fspecial('disk',51),'replicate');


        frameSKEL = bwmorph(frameMSK,'skeleton',inf);
        [skelP(:,1),skelP(:,2)] = find(frameSKEL);

        bp = bwmorph(frameSKEL,'branchpoints');
        [b(:,1),b(:,2)] = find(bp);

        figure(h1);

        BUMP = 11;
        innerMSK = imdilate(innerMSK,strel('square',BUMP));
        paperMSK_s = imerode(paperMSK,strel('square',BUMP));

        paperE = paperMSK_s - innerMSK;
        paperEC = bsxfun(@times,I,paperE);

        BLANK_Paper = I;
        fidxE = find(paperE==1);
        fidxC = find(innerMSK==1);
        for k = 1:3
            tmp = paperEC(:,:,k);
            recolor(k) = mean(tmp(fidxE));
            tmp = BLANK_Paper(:,:,k);
            tmp(fidxC) = recolor(k);
            BLANK_Paper(:,:,k) = tmp;
        end

        HSVn = rgb2hsv(BLANK_Paper);
        HSVn(:,:,3) = GRAD;
        BLANK_Paper = hsv2rgb(HSVn);
        %imshow(BLANK_Paper,[]);
        %waitforbuttonpress




        out = flattenMaskOverlay(I,frameMSK,.6,'r');
        out = flattenMaskOverlay(out,blueSquareMSK,.4,'b');
        out = flattenMaskOverlay(out,frameRegions1,.1,'g');
        out = flattenMaskOverlay(out,frameRegions2,.4,'g');
        out = flattenMaskOverlay(out,frameRegions3,.7,'g');
        %out = flattenMaskOverlay(out,paperMSK_s,.2,'g');




        BLANK_paper = I;



        pt = markFrame(frameMSK,BL,out);
        hold on
        bluePT = [];
        for reg = 1:numel(blueSQUARE_REG)
            de = zeros(4,2);
            ptO = blueSQUARE_REG(reg).BoundingBox(1:2);
            de(1,:) = blueSQUARE_REG(reg).BoundingBox(3:4);
            de(1,1) = 0;
            de(2,:) = blueSQUARE_REG(reg).BoundingBox(3:4);
            de(2,2) = 0;
            de(3,:) = blueSQUARE_REG(reg).BoundingBox(3:4);
            de(4,:) = [0 0];
            de = bsxfun(@plus,de,ptO);
            bluePT(:,:,reg) = de;
            plot(de(:,1),de(:,2),'b*');
        end
        hold off
        drawnow
        save(QRFileList{e},'pt','bluePT','-append');
    catch ME
        ME
    end
    
end
%% work zone
nVCS = bsxfun(@minus,VCS,mean(VCS,2));
[vS vC vU vE vL vERR vLAM] = PCA_FIT_FULL(nVCS,4);
%% mark conetainers MODERN
% might need to go back and re crop with new width and height for conetaners
conetainerData_oPath = '/mnt/tetra/nate/phunSeedlings2/conetainer_EX/';
conetainerData_oPathSIM = '/mnt/tetra/nate/phunSeedlings2/SIMconetainer/';
mkdir(conetainerData_oPathSIM)
conetainerData_oPathSIMW = '/mnt/tetra/nate/phunSeedlings2/SIMconetainer/withCONE/';
mkdir(conetainerData_oPathSIMW)
conetainerData_oPathSIMWO = '/mnt/tetra/nate/phunSeedlings2/SIMconetainer/withoutCONE/';
mkdir(conetainerData_oPathSIMWO)
FilePath = conetainerData_oPath;
coneFileList = {};
FileExt = {'mat'};
coneFileList = gdig(FilePath,coneFileList,FileExt,1);
%sel = [];
FSZ = 101;
VCS = [];
for e = 1:numel(coneFileList)
    try
        sourceFile = coneFileList{e};
        targetFile = strrep(sourceFile,conetainerData_oPath,conetainerData_oPathSIMW);
        copyfile(sourceFile,targetFile);
        targetFile2 = strrep(sourceFile,conetainerData_oPath,conetainerData_oPathSIMWO);
        copyfile(sourceFile,targetFile2);
        
        
        d = load(coneFileList{e});
        I = d.coneImage;
        MASK = d.coneMask;
        size(I,2);
        
        
        msk = d.coneMask;
        msk = imfill(msk,'holes');
        BW = edge(msk);
        BW = imdilate(BW,strel('disk',3,0));

        
        % find the top of the container
        TH = linspace(-15,15,100);
        [H,T,R] = hough(BW','Theta',TH);
        P  = houghpeaks(H,1);
        lines = houghlines(BW',T,R,P,'FillGap',1000,'MinLength',.5*size(I,2));
        k= 1;
        xy = [lines(k).point1; lines(k).point2];

        % calc the rotation parameter and rotate the image and mask
        slope = diff(xy,1,1);
        angle = atan2(slope(1),slope(2));
        I = imrotate(I,angle*180/pi,'bilinear','crop');
        msk = imrotate(msk,angle*180/pi,'bilinear','crop');


        %%%%%%%%%%%%% post rotation
        % find the top of the container
        msk = imfill(msk,'holes');
        BW = edge(msk);
        BW = imdilate(BW,strel('disk',3,0));

        TH = linspace(-15,15,100);
        [H,T,R] = hough(BW','Theta',TH);
        P  = houghpeaks(H,1);
        lines = houghlines(BW',T,R,P,'FillGap',1000,'MinLength',.5*size(I,2));
        k = 1;
        xy = [lines(k).point1; lines(k).point2];
        LEVELLOC  = mean(xy(:,1));
        TALL = round(LEVELLOC);
        

        G = rgb2gray(I);
        E = edge(G);
        sE = sum(E,2);
        sig = mean(G,2);
        figure;
        plot(sig)
        
        
        simfM = imreconstruct(-sig - .05,-sig);
        nsig = -sig - simM;
        TRO = nsig > .05/2;
        TRO(1:20) = 0;
        RT = regionprops(logical(TRO),'Centroid');
        
        
        
        simMU = imreconstruct(sig - .05,sig);
        nsigU = sig - simMU;
        TROU = nsigU > .05/2;
        TROU(1:TALL) = 0;
        RTU = regionprops(logical(TROU),'Centroid');
        
        
        trL = {};
        NT = 4;
        for tp = 1:NT
            trL{tp} = [[1 size(G,2)];[RT(tp).Centroid(2) RT(tp).Centroid(2)]];
        end
        
        trLU = {};
        NTU = 4;
        for tp = 1:NTU
            trLU{tp} = [[1 size(G,2)];[RTU(tp).Centroid(2) RTU(tp).Centroid(2)]];
        end
        
        
        hold on
        plot(nsig,'b')
        hold off
        close all
        %waitforbuttonpress
        
        
        
        
        sig = im2col(sig,[FSZ 1],'sliding');


        for s = 1:size(sig,2)
            sig(:,s) = sig(:,s) / norm(sig(:,s));
        end

        if e < -50
            [l1,l2,~] = impixel(G);
            l2 = l2 -(FSZ-1)/2;
            sel = [sel sig(:,l2)];
            
        end
        
        sigBLOCK = sig;
        
        
        sig = mean(sel,2)'*sig;
        sig(round(xy(:,2) + 200):end) = 0;
        [~,sidx] = max(sig);
        
        %sel = [sel sigBLOCK(:,sidx)];
        
        sidx
        l2
        sidx = sidx + (FSZ-1)/2 - 30;

        
        
        %{
        if sidx > 50
            msk(sidx-50:end,:) = 1;
        end
        %}
        
        
       


       
        %cLOC = find(mean(msk(1:TALL,:),1) > .35);
        mskL = imopen(msk,strel('disk',21,0));
        cLOC = find(mean(mskL(1:TALL,:),1) > .6*max(mean(mskL(1:TALL,:),1)));
        %{
        if LEVELLOC > 9
            SNIP = [9 5];
            cLOC = find(mean(msk(LEVELLOC-SNIP(1):LEVELLOC-SNIP(2),:),1) > .2);
        end
        %}
        [gG U8] = gradient(G);
        gG = abs(gG);
        SG = std(G,1,2);
        figure;
        sig = std(gG,1,2);
      
        sU = mean(abs(U8),2);
        sig = sU.*sig;
        
        %sig = gradient(SG);
        kS = 11;
        sig = imfilter(sig,ones(kS,1)/kS,'replicate');
        
       
        LEVEL = .1*10^-4;
        simMU = imreconstruct(sig - LEVEL,sig);
        nsigU = sig - simMU;
        TROU = nsigU > LEVEL/2;
        TROU(1:(trLU{end}(2,1)+21)) = 0;
        RTU2 = regionprops(logical(TROU),'Centroid');
        
        %figure;plot(TROU);
        %waitforbuttonpress
        
        trLU2 = {};
        NTU2 = 3;
        for tp = 1:NTU2
            trLU2{tp} = [[1 size(G,2)];[RTU2(tp).Centroid(2) RTU2(tp).Centroid(2)]];
        end
       
        

        imshow(G,[]);
        hold on
        
        
        for tp = 1:numel(trL)
            plot(trL{tp}(1,:),trL{tp}(2,:),'y')
        end
        
        
        for tp = 1:numel(trLU)
            plot(trLU{tp}(1,:),trLU{tp}(2,:),'c')
            FY = mean(trLU{tp}(2,:));
        end
        CONE_WID = 30;
        for tp = 1:numel(trLU2)
            plot(trLU2{tp}(1,:),trLU2{tp}(2,:),'m');
            if tp == 2
                CONESTR = G(trLU2{1}(2,1):trLU2{3}(2,1),:);
                CONESTR = imresize(CONESTR,[100 1500]);
            end
        end
        %{
        nCONESTR = bsxfun(@minus,CONESTR',mean(CONESTR',2));
        vC = PCA_REPROJ(nCONESTR,vE,vU);
        SIMY = PCA_BKPROJ(vC,vE,vU);
        SIMY = bsxfun(@plus,SIMY,mean(CONESTR',2));
        close all
        plot(vC);
        waitforbuttonpress
        VCS = [VCS;CONESTR'];
        %}
        
        
        %{
        imshow(CONESTR,[]);
        
        
        [cg1,cg2] = gradient(CONESTR);
        cg2 = imfilter(cg2,fspecial('gaussian',[11 11],5),'replicate');
        [ccg1,ccg2] = gradient((cg2));
        GGG = (ccg1);
        plot(mean(GGG,1))
        
        
        cg1 = abs(cg1);
        imshow(cg2,[]);
        waitforbuttonpress
        %}
        IY = mean(xy(:,1));
        
        CPPP = xy(1,:) + diff(xy,1,1)/2;
        
        
        DY = (FY - IY)/10;
        DX = linspace(-130,130,10);
        SLOPE = (DX).*DY.^-1;
        [g1,g2] = ndgrid(linspace(0,FY-IY,10),linspace(-FY-IY/2,FY-IY/2,10));
        
        SLOPE = repmat(SLOPE,[size(g1,1) 1]);
        XP = g1.*SLOPE;
        plot(XP(:)+CPPP(2),g1(:)+CPPP(1))
        
        
        
        plot(1:size(G,2),sidx*ones(1,size(G,2)),'r');
        plot(cLOC(1)*ones(1,size(G,1)),1:size(G,1),'r');
        plot(cLOC(end)*ones(1,size(G,1)),1:size(G,1),'r');
        
        
        plot(CPPP(2),CPPP(1),'m^')
        
        
        
        
        points = detectHarrisFeatures(G);
        %plot(points.Location(:,1),points.Location(:,2),'r*');
        xy(:,2) = [1;size(G,2)]';
        plot(xy(:,2),xy(:,1),'g');

        
        newMask = MASK;
        newMask(1:xy(:,1),:) = 0;
        newMask(round(trLU2{end}(2,2)):end,cLOC(1):cLOC(end)) = 0;
        %imshow(newMask,[])
        
        
        hold off
        drawnow
        coneMask = newMask;
        save(targetFile2,'coneMask','-append');
        %waitforbuttonpress
        pause(.3)
        WIDTH(e) = cLOC(end) - cLOC(1);
        HEIGHT(e) = mean(points.Location(:,2)) - sidx;
    catch ME
        ME
        
    end
end
%% load marks for QR

qrcode_oPath = '/mnt/tetra/nate/phunSeedlings/qrcodes/';
FilePath = qrcode_oPath;
QRFileList = {};
FileExt = {'mat'};
QRFileList = gdig(FilePath,QRFileList,FileExt,1);

qrNORSZ = [796 966];
scale = .25;
qrSZ = round(scale*[796 966]);

patchSZ = round(scale*[51 51]);
e =1;
numToLoad = 750;
d = load(QRFileList{e});
I = double(d.QR);
oSZ = size(I);
I = imresize(I,qrNORSZ);
toScale = oSZ.*size(I).^-1;
I = imresize(I,qrSZ);


maskSZ = size(I);
maskSZ(end) = [];
pointList = round(d.pt);

da = 3;
[tmpMask] = generatePointMasks(maskSZ,pointList,scale,da);
PERR = .005;
[tmpMaskBLUE] = generateBlueBoxMarks(maskSZ,round(d.bluePT),scale*toScale(1:2),da);
tmpMask = cat(3,tmpMask,tmpMaskBLUE);
[tmp_squareX,tmp_squareY] = generateSquarescanData(I,tmpMask,patchSZ,[PERR,1]);
%[tmp_squareXB,tmp_squareYB] = generateSquarescanData(I,tmpMaskBLUE,patchSZ,[PERR,1]);

szX = size(tmp_squareX);
szY = size(tmp_squareY);
szX(4) = szX(4)*numToLoad;
szY(1) = szY(1)*numToLoad;
qrX = zeros(szX);
qrY = zeros(szY);


%szXB = size(tmp_squareXB);
%szYB = size(tmp_squareYB);
%szXB(4) = szXB(4)*numToLoad;
%szYB(1) = szYB(1)*numToLoad;
%qrXB = zeros(szXB);
%qrYB = zeros(szYB);

SKIP = size(tmp_squareX,4);

%SKIPB = size(tmp_squareXB,4);
str = 1;
strB = 1;
loadCNT = 1;
ptr = 1;

while loadCNT <= numToLoad
    try
        
       
        d = load(QRFileList{ptr});
        ptr = ptr + 1;
        I = double(d.QR);
        oSZ = size(I);
        I = imresize(I,qrNORSZ);
        toScale = size(I).*oSZ.^-1;
        I = imresize(I,qrSZ);
        
        maskSZ = size(I);
        maskSZ(end) = [];
        pointList = round(d.pt);
        if ~any(isnan(pointList(:)))
            stp = str + SKIP - 1;
            %stpB = strB + SKIPB - 1;
            
            
            [tmpMask] = generatePointMasks(maskSZ,pointList,scale*toScale(1:2),da);
            [tmpMaskBLUE] = generateBlueBoxMarks(maskSZ,round(d.bluePT),scale*toScale(1:2),da);
            tmpMask = cat(3,tmpMask,tmpMaskBLUE);
            
            [tmp_squareX,tmp_squareY] = generateSquarescanData(I,tmpMask,patchSZ,[PERR,1]);
            %[tmp_squareXB,tmp_squareYB] = generateSquarescanData(I,tmpMaskBLUE,patchSZ,[PERR,1]);
            
            
            qrX(:,:,:,str:stp) = tmp_squareX;
            qrY(str:stp,:) = tmp_squareY;
            %qrXB(:,:,:,strB:stpB) = tmp_squareXB;
            %qrYB(strB:stpB,:) = tmp_squareYB;
            
            
            
            str = stp + 1;
            %strB = stpB + 1;
            
            loadCNT = loadCNT + 1;
        end
        
    catch MEa
        
        MEa
    end
    loadCNT
end
%%
zidx = find(all(qrY==0,2));
zidx1 = find(any(qrY==1,2));
extraC = zeros(size(qrY,1),1);
extraC(zidx) = 1;
NqrY = [qrY extraC];
sum(NqrY,2);
NqrY = vec2ind(NqrY');
%% train QR location regression
imgSZ = size(qrX);
imgSZ(4) = [];
inputLayer = imageInputLayer([imgSZ]);

middleLayer1 = [...
    convolution2dLayer([7 7],15)
    reluLayer
    crossChannelNormalizationLayer(4)
    maxPooling2dLayer([2 2],'Stride',1)];
middleLayer2 = [...
    convolution2dLayer([3 3],8)
    reluLayer
    crossChannelNormalizationLayer(4)
    maxPooling2dLayer(2,'Stride',2)];
finalLayers = [...
    fullyConnectedLayer(9+4)
    softmaxLayer
    classificationLayer];
layers_mark = [...
    inputLayer
    middleLayer1
    middleLayer2
    finalLayers];
options = trainingOptions('sgdm',...
    'MaxEpochs',1000, ...
    'Verbose',true,...
    'InitialLearnRate',.005,...
    'Plots','training-progress',...
    'L2Regularization',.0001,...
    'MiniBatchSize',128*4,...
    'ExecutionEnvironment','parallel');
%% training count network
trainedNet_QR_mark2 = trainNetwork(qrX,categorical(NqrY'),layers_mark,options);
%% try testing QR marks
for LL = 1:1%numel(QRFileList)
    d = load(QRFileList{LL});
    I = double(d.QR);
    oSZ = size(I);
    I = imresize(I,qrNORSZ);
    toScale = size(I).*oSZ.^-1;
    oSZ2 = size(I);
    I = imresize(I,qrSZ);
    toScale2 = size(I).*oSZ2.^-1;
    [markX,tmp_squareY] = generateSquarescanData(I,[],patchSZ,[]);
    MARKS = trainedNet_QR_mark2.classify(markX,'MiniBatchSize',256*512);
    V = reshape(double(MARKS(:,1)),[size(I,1) size(I,2)]-(patchSZ-1));
    imshow(V,[])
end
%% try testing QR location
qrSZ = round(.15*[796 966]);
qrI_X = [];
qrI_Y = [];
for e = 400:numel(QRFileList)
    try
        d = load(QRFileList{e});
        TQR = imresize(d.QR,qrSZ);
        LOCT = trainedNetLOCR.predict(TQR);
        LOCT = reshape(LOCT,[8 2]);
        imshow(d.QR,[]);
        hold on
        plot(LOCT(:,2),LOCT(:,1),'k*')
        hold off
        drawnow
    catch ME
        ME
    end
end
%% load QR data for clustering
clear QR
M = [];
for e = 1:5:numel(QRFileList)
    load(QRFileList{e},'QR')
    QR = imresize(QR,.25);
    sz = size(QR);
    QR = reshape(QR,[prod(sz(1:2)) sz(3)]);
    M = [M;QR];
    e
end
%% for labeling QR images
QRlabel = gmdistribution.fit(M,4);
QRlabel2_2 = gmdistribution.fit(M,5);
QRlabel3 = gmdistribution.fit(M,6);
%% view cluster
close all
for e = 1:30
 
    load(FileList{e},'QR')
    sz = size(QR);
    I = QR;
    GI = rgb2gray(I);
    Cpoints = detectHarrisFeatures(GI);
    
    QR = reshape(QR,[prod(sz(1:2)) sz(3)]);
    k = QRlabel2.cluster(QR);
    k = reshape(k,sz(1:2));
    rgb = label2rgb(k);
    imshow(rgb,[]);
    drawnow
end
%% model QR sheet
close all
clear toConnect
sel = [2 1];
sel = [5 3 2];
sel = [1 3 5];
sel = [1 3 5];
h1 = figure;
h2 = figure;
boxCNT = 1;
blueBOX_stack = [];
blueBOX_size = [];
LOC = [];
frameY = [];
frameX = [];
for e = 1:numel(FileList)
    
    
    clear b PO skelP
    load(FileList{e},'QR')
    sz = size(QR);
    I = QR;
    GI = rgb2gray(I);
    Cpoints = detectHarrisFeatures(GI);
    
    QR = reshape(QR,[prod(sz(1:2)) sz(3)]);
    [k,~,KK] = QRlabel2_2.cluster(QR);
    k = reshape(k,sz(1:2));
    KK = reshape(KK,[sz(1:2) size(KK,2)]);
    rgb = label2rgb(k);
    frameMSK = k == sel(1);
    paperMSK = k ~= sel(2);
    blueSquareMSK = KK(:,:,sel(3)) > .98;
    blueSquareMSK = bwareaopen(blueSquareMSK,100);
   
    %blueSquareMSK = bwlarge(blueSquareMSK,24);
    frameMSK = bwareaopen(frameMSK,300);
    frameMSK = imclose(frameMSK,strel('line',20,0));
    frameMSK = imclose(frameMSK,strel('line',20,90));
    %frameMSK = imerode(frameMSK,strel('square',5));
    frameMSK = bwlarge(frameMSK);
    
    
    frameRegions = imfill(frameMSK,'holes') - frameMSK;
    frameRegions1 = logical(bwlarge(frameRegions,1));
    frameRegions2 = logical(bwlarge(frameRegions,2) - frameRegions1);
    frameRegions3 = logical(bwlarge(frameRegions,3) - frameRegions2 - frameRegions1);
    
    
    
    blueSquareMSK = logical(blueSquareMSK.*frameRegions1);
    
    blueSquareMSK = imfill(blueSquareMSK,'holes');
    blueSquareMSK = bwlarge(blueSquareMSK,24);
    blueSQUARE_REG = regionprops(blueSquareMSK,'PixelIdxList','BoundingBox','Centroid');
    
    BOXCEN = [];
    for box = 1:numel(blueSQUARE_REG)
        BOXCEN(box,:) = blueSQUARE_REG(box).Centroid;
        img = imcrop(I,blueSQUARE_REG(box).BoundingBox);
        blueBOX_size(boxCNT,:) = [size(img,1) size(img,2)];
        img = imresize(img,[65 80]);
        blueBOX_stack(:,:,:,boxCNT) = img;
        figure(h2);
        imshow(img,[]);
        boxCNT = boxCNT + 1;
        drawnow
    end
    
    nBOXSTACK = [];
    for col = 1:6
        [~,bidx] = sort(BOXCEN(:,1));
        [~,bidx2] = sort(BOXCEN(bidx(1:4),2));
        nBOXSTACK = [nBOXSTACK;BOXCEN(bidx(bidx2),:)];
        BOXCEN(bidx(bidx2),:) = [];
    end
    BOXCEN = nBOXSTACK;
    
    
    paperMSK = bwlarge(paperMSK);
    paperMSK = imfill(paperMSK,'holes');
    paperR = regionprops(paperMSK);
    
    
    innerMSK = imfill(frameMSK,'holes');
    
    
    LIGHTING = k == 4;
   
    HSV = rgb2hsv(I);
    GRAD = imfilter(HSV(:,:,3),fspecial('disk',51),'replicate');
    
    
    frameSKEL = bwmorph(frameMSK,'skeleton',inf);
    [skelP(:,1),skelP(:,2)] = find(frameSKEL);
    
    bp = bwmorph(frameSKEL,'branchpoints');
    [b(:,1),b(:,2)] = find(bp);
    
    figure(h1);
    
    BUMP = 11;
    innerMSK = imdilate(innerMSK,strel('square',BUMP));
    paperMSK_s = imerode(paperMSK,strel('square',BUMP));
    
    paperE = paperMSK_s - innerMSK;
    paperEC = bsxfun(@times,I,paperE);
    
    BLANK_Paper = I;
    fidxE = find(paperE==1);
    fidxC = find(innerMSK==1);
    for k = 1:3
        tmp = paperEC(:,:,k);
        recolor(k) = mean(tmp(fidxE));
        tmp = BLANK_Paper(:,:,k);
        tmp(fidxC) = recolor(k);
        BLANK_Paper(:,:,k) = tmp;
    end
   
    HSVn = rgb2hsv(BLANK_Paper);
    HSVn(:,:,3) = GRAD;
    BLANK_Paper = hsv2rgb(HSVn);
    %imshow(BLANK_Paper,[]);
    %waitforbuttonpress
    
    
    
    
    out = flattenMaskOverlay(I,frameMSK,.6,'r');
    out = flattenMaskOverlay(out,blueSquareMSK,.4,'b');
    out = flattenMaskOverlay(out,frameRegions1,.1,'g');
    out = flattenMaskOverlay(out,frameRegions2,.4,'g');
    out = flattenMaskOverlay(out,frameRegions3,.7,'g');
    %out = flattenMaskOverlay(out,paperMSK_s,.2,'g');
    
    
    
    
    BLANK_paper = I;
    
    
    
    pt = markFrame(frameMSK,BL,out);
    hold on
    for p = 1:size(pt,1)
        [idx(p)] = snapTo(skelP,pt(p,:));
        PO(p,:) = skelP(idx(p),:);
    end
    
    for p = 1:size(PO,1)
        plot(PO(p,2),PO(p,1),'ko');
        text(PO(p,2),PO(p,1),num2str(p),'BackGroundColor','y');
    end
    
    CP = fliplr(PO(7,:));
    plot(CP(1),CP(2),'b*');
    %waitforbuttonpress
    tmpLOC = [bsxfun(@minus,fliplr(PO),CP);bsxfun(@minus,BOXCEN,CP)];
    LOC(:,:,e) = tmpLOC;
    tmpLOCsz = size(tmpLOC);
    tmpLOC = reshape(tmpLOC,[prod(tmpLOCsz(1:2)) 1]);
    
    lC = PCA_REPROJ_T(tmpLOC,lE,lU);
    lSIM = PCA_BKPROJ_T(lC,lE,lU);
    lSIM = reshape(lSIM,tmpLOCsz);
    tmpLOC = reshape(tmpLOC,tmpLOCsz);
    lSIM = bsxfun(@plus,lSIM,CP);
    tmpLOC = bsxfun(@plus,tmpLOC,CP);
    
    plot(lSIM(:,1),lSIM(:,2),'r.');
    plot(tmpLOC(:,1),tmpLOC(:,2),'b.');
    
    RR = Radjacency(skelP', 2);
    
    toConnect{1} = [[3 5];[5 7];[7 8];[8 3]];
    toConnect{2} = [[5 1];[1 6];[6 7];[7 5]];
    toConnect{3} = [[8 6];[6 2];[2 4];[4 8]];
    iLEN{1} = 100;
    iLEN{2} = 300;
    iLEN{3} = 500;
    TOTP = [];
    TMPframeY = [];
    TMPframeX = [];
    for r = 1:numel(toConnect)
        for t = 1:size(toConnect{r},1)
            [path , pathcost] = dijkstra(RR , idx(toConnect{r}(t,1)) , idx(toConnect{r}(t,2)));
            iPATH = interp1(1:numel(path),skelP(path,:),linspace(1,numel(path),iLEN{r}));
            
            [DomainS,DomainG] = extendChromosomeMidline(iPATH,[0 0],I,2);

            %size(iPATH)
            %waitforbuttonpress
            dsz = size(DomainG);
            vec = [];
            for k = 1:size(I,3)
                vec(:,k) = ba_interp2(I(:,:,k),DomainS(:,1),DomainS(:,2));
            end
            vec = reshape(vec,[dsz(1) dsz(2) 3]);
            
            FX = size(DomainG);
            
            DomainS = bsxfun(@minus,DomainS,CP);
            
            TMPframeY = [TMPframeY;vec(:)];
            TMPframeX = [TMPframeX;DomainS(:)];
            
            
            %{
            for s1 = 1:5:size(DomainG,1)
                plot(DomainG(s1,:,1),DomainG(s1,:,2),'k')
            end
            for s1 = 1:5:size(DomainG,2)
                plot(squeeze(DomainG(:,s1,1)),squeeze(DomainG(:,s1,2)),'k')
            end
            drawnow
            %}
            
            
            TOTP = [TOTP;iPATH];
            plot(skelP(path,2),skelP(path,1),'k')
        end
    end
    FX = size(TMPframeX);
    frameX = [frameX TMPframeX(:)];
    frameY = [frameY TMPframeY(:)];
    
    innerW = skelP(idx(3),2) - skelP(idx(1),2);
    innerH = skelP(idx(3),1) - skelP(idx(4),1);
    
    innerBOX = [fliplr(skelP(idx(3),:)) -innerW -innerH];
    EX_x = (sum(any(paperMSK,1)) + innerW)/4;
    EX_y = (sum(any(paperMSK,2)) + innerH)/4;
    PER = .8;
    innerBOX(1:2) = innerBOX(1:2) - PER*[EX_x EX_y];
    innerBOX(3:4) = innerBOX(3:4) + PER*2*[EX_x EX_y];
    
    
    
    
    %plot(b(:,2),b(:,1),'b*');
    %rectangle('Position',paperR.BoundingBox)
    %plot(Cpoints.Location(:,1),Cpoints.Location(:,2),'ko')
    %rectangle('Position',innerBOX);
    
    %waitforbuttonpress
    hold off
    figure(h2);
    imshow(rgb,[]);
    drawnow
    %waitforbuttonpress
    e
end
%% decompose the blue squares
close all
boxSZ = size(blueBOX_stack);
data = reshape(blueBOX_stack,[prod(boxSZ(1:3)) boxSZ(4)]);
[qU,qE,qL] = PCA_FIT_FULL_Tws(data,3);
qC = PCA_REPROJ_T(data,qE,qU);
plot3(qC(1,:),qC(2,:),qC(3,:),'.');
qSIM = PCA_BKPROJ_T(qC,qE,qU);
qSIM = reshape(qSIM,boxSZ);
%% decompose loc data
close all
locSZ = size(LOC);
data = reshape(LOC,[prod(locSZ(1:2)) locSZ(3)]);
[lU,lE,lL] = PCA_FIT_FULL_Tws(data,3);
lC = PCA_REPROJ_T(data,lE,lU);
plot3(lC(1,:),lC(2,:),lC(3,:),'.');
lSIM = PCA_BKPROJ_T(lC,lE,lU);
lSIM = reshape(lSIM,locSZ);
%%
[BL] = makeFramePipe();
%% try it on whole images and gather container-black and shadow blue and plant and white from container
close all
cnt = 1;
HV = [];
toAdd = [];
totD = 200;

for fl = 1:numel(imageFileList)
    for e = 1:totD
        try
            %if strcmp(imageFileList{fl}{e}(end-4),'}')
            if ~contains(imageFileList{fl}{e},'output')
                
                I = double(imread(imageFileList{fl}{e}))/255;

                [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(I,20);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% find background
                % look at the size of the crop box
                sz = size(bioStrip);
                % reshape the data
                tmp1 = reshape(bioStrip,[prod(sz(1:2)) sz(3)]);

                [K,~,nl] = obj.cluster(tmp1);
                [K] = reshape(K,sz(1:2));
                [nl] = reshape(nl,[sz(1:2) size(nl,2)]);
                BK = K == 1;
                %%%%%%%% find background
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
               
                bkTMP = bsxfun(@times,bioStrip,BK);
                S = (sum(bkTMP,1));
                M = (sum(BK,1));
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% find foreground
                FG = K == 2;
                newPot = find(FG);
                newPot = newPot(randperm(numel(newPot)));
                W = [];
                for k = 1:size(bioStrip,3)
                    tmp = bioStrip(:,:,k);
                    W = [W tmp(newPot(1:800))];
                end
                toAdd = [toAdd;W];
                %%%%%%%% find foreground
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %{
                newPot = find(nl(:,:,2) < .95 & K ~= 1);
                newPot = newPot(randperm(numel(newPot)));
                W = [];
                for k = 1:size(bioStrip,3)
                    tmp = bioStrip(:,:,k);
                    W = [W tmp(newPot(1:200))];
                end
                toAdd = [toAdd;W];
                %}
                
                X = linspace(1,size(I,2),5000);
                M = interp1(1:size(I,2),M,X);
                for k = 1:size(S,3)
                   nS(:,:,k) = interp1(1:size(I,2),S(:,:,k),X);
                end
                S = nS;
                HV(cnt,:,:) = bsxfun(@times,S,M.^-1);
                
                cnt = cnt + 1;
                
                out = flattenMaskOverlay(bioStrip,bwlarge(K==1),.5,'r');
                imshow(out,[]);
                drawnow
                %waitforbuttonpress
            end
        catch ME
            ME
            %break
        end
        e
    end
    
end
%% cluster plants from background and container
nextLevelCluster = [toAdd;colorGroups{2}];
nonBackgroundGMM = fitgmdist(nextLevelCluster,5);
%% try sub clustering
coneTainerModelDepth = 400;
totD = 300;

coneTainerStack = zeros(coneTainerModelDepth,1500,3,totD*numel(imageFileList));
coneTainerStackPOS = zeros(totD*numel(imageFileList),1);
coneTainerLocation = zeros(totD*numel(imageFileList),2);
coneTainerWidth = zeros(totD*numel(imageFileList),1);

qrStack = zeros(600,775,3,totD*numel(imageFileList));

plantBufferWidth = 125;
plantStackColor = zeros(2000,1500,3,totD*numel(imageFileList));
plantStackMask = zeros(2000,1500,totD*numel(imageFileList));
cnt = 1;
Qcnt = 1;
oPath = '/mnt/tetra/nate/sceneModel/';
oPath = '';
totD = 400;


for fl = 1:numel(imageFileList)
    parfor e = 1:totD
        try
            %if strcmp(imageFileList{fl}{e}(end-4),'}')
            if ~contains(imageFileList{fl}{e},'output')
                
                generateSceneData(imageFileList{fl}{e},obj,nonBackgroundGMM,[3000 1500],[400 1500],[600 775],125,20,400,oPath);
   
                
                %waitforbuttonpress
            end
        catch ME
            ME
            %break
        end
        e
    end
    
end
%% try cluster and sub cluster on crop boxes
for fl = 1:numel(FileList)
    for e = 1:totD
        try
            tmp = load(FileList{fl}{e});
            % for each plant
            for i = 1:numel(tmp.returnI)
                
                % look at the size of the crop box
                sz = size(tmp.returnI{i});
                % reshape the data
                tmp1 = reshape(tmp.returnI{i},[prod(sz(1:2)) sz(3)]);
                
                K = obj.cluster(tmp1);
                K = reshape(K,sz(1:2));
                
                out = flattenMaskOverlay(tmp.returnI{i},bwlarge(K==2),.5,'r');
                out = flattenMaskOverlay(out,tmp.MASK{i}==1,.5,'b');
                imshow(out,[]);
                drawnow
              

            end
        catch ME
            ME
            %break
        end
        e
    end
    
end
%% scan for tables and ingest them
FilePath = '/mnt/tetra/nate/seedlingImageParts/dataTables/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
masterTable = table;
for e = 1:numel(FileList)
    %tmp = readtable(FileList{e});
    tmp = load(FileList{e});
    masterTable = [masterTable;tmp.opTable];
end
%%
modelTable = table;
modelTable = generateRegConetainerData(masterTable,modelTable);
%% label the QR
IDX = findImageFramePair(masterTable,'qr_image');
newMeta = generateRegQRData(masterTable(IDX,:),'/mnt/tetra/nate/modelData/GMM.mat');
masterTable = [masterTable;newMeta];
%%
rawIDX = find(strcmp(masterTable.type,'maizeSeedling_rawImage'));
rawKEY = masterTable.key(rawIDX);
[IDX] = getAllGeometryObjectsForImageKey(masterTable,rawKEY(1));

%%
close all
imshow(simImage,[])
hold all
POS = generateObjectDistributions(masterTable);
for e = 1:numel(POS.CM)
    plot(POS.CM{e}(:,2),POS.CM{e}(:,1),'.')
    hold all
    LABEL{e} = strrep(POS.LABEL{e},'_','');
end
legend(LABEL)
%% break some into multiple distributions
toSearch = 'qrObject';
fidx = find(strcmp(POS.LABEL,toSearch));
k = 1;
qrDistribution = fitgmdist(POS.LOC{fidx},k);

toSearch = 'conetainer_whole';
fidx = find(strcmp(POS.LABEL,toSearch));
k = 3;
containerLocation = fitgmdist(POS.LOC{fidx},k);
%%
for e = 1:100
    sSZ = [3280 4992]; 
    simImage = generateScene(masterTable,qrDistribution,containerLocation,sSZ);
    close all
    imshow(simImage,[]);
    drawnow
end
%% model QR
qrFilePath{1} = '/mnt/tetra/nate/seedlingImageParts/qrCodes/';
for fl = 1:numel(qrFilePath)
    qrFileList{fl} = {};
    FileExt = {'tif'};
    qrFileList{fl} = gdig(qrFilePath{fl},qrFileList{fl},FileExt,1);
end
%% make pipe blocks
close all
W = 15;
w = 4;
B = zeros(2*W+1);
B((end-1)/2+1,(end-1)/2+1) = 1;
B((end-1)/2-w+1:(end-1)/2+w+1,1:((end-1)/2+1+w)) = 1;
imshow(B,[]);

BB = [];
BB(:,:,1) = B;
for e = 2:4
    BB(:,:,e) = imrotate(BB(:,:,e-1),90);
    imshow(BB(:,:,e),[]);
    drawnow
end
BL = [];
cnt = 1;
for e = 2:4
    v = nchoosek(1:4,e);
   v
    for k = 1:size(v,1)
        newB = zeros(size(B));
        for i = 1:size(v,2)
            newB = newB | (BB(:,:,v(k,i)) == 1);
        end
        newB = newB .* sum(newB(:)).^-1 + e*.00001;
        BL(:,:,cnt) = newB;
        
        imshow(BL(:,:,cnt),[])
        drawnow
        %waitforbuttonpress
        cnt = cnt + 1;
       
    end
    
end

%%
close all
boxSZ = 6;

I1 = imread(qrFileList{1}{2});
newM = [];
newI = [];
registeredQRpath = '/mnt/tetra/nate/seedlingImageParts/qrCodes/regQRcode/';

for i = 1:numel(qrFileList{1})
    [~,nm] = fileparts(qrFileList{1}{i});
    I = imread(qrFileList{1}{i});
    sz = size(I);
    CI = reshape(I,[prod(sz(1:2)) sz(3)]);
    T = GMModel.cluster(double(CI));
    T = reshape(T,sz(1:2));
    %{
    G1 = rgb2gray(I1);
    points = detectMinEigenFeatures(G1);
    %}
    T = bwlarge(T==2);
    T = imclose(T,strel('square',7));
    %T = imfill(T,'holes');
    R = [];
    for c = 1:size(BL,3)
        R(:,:,c) = imfilter(double(T),double(BL(:,:,c)),'replicate');
    end
    [IDXP,IDX] = max(R,[],3);
    eT = imerode(T,strel('square',5));
    mIDX = IDX.*eT;
    %imshow(mIDX,[]);
    toSearchFor = [1 3 4 6 7 8 9 10];
    pt = [];
    for e = 1:numel(toSearchFor)
        ff = bwlarge(mIDX == toSearchFor(e));
        fidx = find(ff);
        [mm,midx] = max(IDXP(fidx));
        midx = find(IDXP(fidx) == mm);
        AL = [];
        [AL(:,1) AL(:,2)] = ind2sub(size(T),fidx(midx));
        pt(e,:) = mean(AL,1);
    end
    
    
    
    imshow(I,[]);
    hold on
    plot(pt(:,2),pt(:,1),'k*')
    drawnow

    
    hold on
    for e = 1:size(pt,1)
        plot(pt(e,2),pt(e,1),'k*');
        text(pt(e,2)+3,pt(e,1),num2str(e));
    end
    
    checkBOX = pt([8 6 4 2],:);
    qrBOX = pt([3 5 8 7],:);
    humanBOX = pt([5 1 7 6],:);
    plot(checkBOX(1:2,2),checkBOX(1:2,1),'k')
    plot(checkBOX(3:4,2),checkBOX(3:4,1),'k')
    plot(checkBOX([1 3],2),checkBOX([1 3],1),'k')
    plot(checkBOX([2 4],2),checkBOX([2 4],1),'k')
    
    plot(qrBOX(1:2,2),qrBOX(1:2,1),'k')
    plot(qrBOX(3:4,2),qrBOX(3:4,1),'k')
    plot(qrBOX([1 3],2),qrBOX([1 3],1),'k')
    plot(qrBOX([2 4],2),qrBOX([2 4],1),'k')
    
    plot(humanBOX(1:2,2),humanBOX(1:2,1),'k')
    plot(humanBOX(3:4,2),humanBOX(3:4,1),'k')
    plot(humanBOX([1 3],2),humanBOX([1 3],1),'k')
    plot(humanBOX([2 4],2),humanBOX([2 4],1),'k')
    drawnow
    pause(.1)
    %waitforbuttonpress
    
    pt = fliplr(pt);
    
    close all
    buildModel = 1;
    if i ~= 1 | buildModel
        
        tform = fitgeotrans(pt,fixedPoints,'similarity');
        %tform = fitgeotrans(fixedPoints,pt,'similarity');
        
        tmpI = imwarp(I,tform,'OutputView',outP);
        newI(:,:,:,i) = tmpI;
        
        newM(:,:,i) = imwarp(T,tform,'OutputView',outP);
        toV = mean(newI/255,4);
        toM = mean(newM,3);
        
        
        out = flattenMaskOverlay(toV,toM > .3,.8,'g');
        
        imshow(out,[]);
        drawnow
        hold off
        %waitforbuttonpress
    else
        fixedPoints = pt;
        outSize = size(I);
        outP = imref2d(outSize(1:2));
        newI(:,:,:,1) = I;
    end
    
    
    if buildModel
        imwrite(newI(:,:,:,i),[registeredQRpath nm '.tif']);
    end
    
    %{
    imshow(I,[]);
    hold on
    plot(pt(:,2),pt(:,1),'k*')
    drawnow

    hold on
    for e = 1:size(pt,1)
        plot(pt(e,2),pt(e,1),'k*');
        text(pt(e,2)+3,pt(e,1),num2str(e));
    end
    
    checkBOX = pt([8 6 4 2],:);
    qrBOX = pt([3 5 8 7],:);
    humanBOX = pt([5 1 7 6],:);
    plot(checkBOX(1:2,2),checkBOX(1:2,1),'k')
    plot(checkBOX(3:4,2),checkBOX(3:4,1),'k')
    plot(checkBOX([1 3],2),checkBOX([1 3],1),'k')
    plot(checkBOX([2 4],2),checkBOX([2 4],1),'k')
    
    plot(qrBOX(1:2,2),qrBOX(1:2,1),'k')
    plot(qrBOX(3:4,2),qrBOX(3:4,1),'k')
    plot(qrBOX([1 3],2),qrBOX([1 3],1),'k')
    plot(qrBOX([2 4],2),qrBOX([2 4],1),'k')
    
    plot(humanBOX(1:2,2),humanBOX(1:2,1),'k')
    plot(humanBOX(3:4,2),humanBOX(3:4,1),'k')
    plot(humanBOX([1 3],2),humanBOX([1 3],1),'k')
    plot(humanBOX([2 4],2),humanBOX([2 4],1),'k')
    %}
end
%%

%% build model

toV = mean(newI/255,4);
toM = mean(newM,3);
frameModel = bsxfun(@times,(toM > .5),toV);
T = (toM > .5);
R = [];
for c = 1:size(BL,3)
    R(:,:,c) = imfilter(double(T),double(BL(:,:,c)),'replicate');
end
[IDXP,IDX] = max(R,[],3);
eT = imerode(T,strel('square',5));
mIDX = IDX.*eT;
%imshow(mIDX,[]);
toSearchFor = [1 3 4 6 7 8 9 10];
pt = [];
for e = 1:numel(toSearchFor)
    ff = bwlarge(mIDX == toSearchFor(e));
    fidx = find(ff);
    [mm,midx] = max(IDXP(fidx));
    midx = find(IDXP(fidx) == mm);
    AL = [];
    [AL(:,1) AL(:,2)] = ind2sub(size(T),fidx(midx));
    pt(e,:) = mean(AL,1);
end

close all
imshow(frameModel,[]);
hold on
for e = 1:size(pt,1)
    plot(pt(e,2),pt(e,1),'k*');
    text(pt(e,2)+3,pt(e,1),num2str(e));
end
checkBOX = pt([8 6 4 2],:);
qrBOX = pt([3 5 8 7],:);
humanBOX = pt([5 1 7 6],:);
[a(1)] = measurePointBox(checkBOX);
[a(2)] = measurePointBox(qrBOX);
[a(3)] = measurePointBox(humanBOX);
toV = imrotate(toV,mean(a)*180/pi);
toM = imrotate(toM,mean(a)*180/pi);

frameModel = bsxfun(@times,(toM > .5),toV);
imshow(toV,[]);
drawnow
T = (toM > .5);
R = [];
for c = 1:size(BL,3)
    R(:,:,c) = imfilter(double(T),double(BL(:,:,c)),'replicate');
end
[IDXP,IDX] = max(R,[],3);
eT = imerode(T,strel('square',5));
mIDX = IDX.*eT;
%imshow(mIDX,[]);
toSearchFor = [1 3 4 6 7 8 9 10];
pt = [];
for e = 1:numel(toSearchFor)
    ff = bwlarge(mIDX == toSearchFor(e));
    fidx = find(ff);
    [mm,midx] = max(IDXP(fidx));
    midx = find(IDXP(fidx) == mm);
    AL = [];
    [AL(:,1) AL(:,2)] = ind2sub(size(T),fidx(midx));
    pt(e,:) = mean(AL,1);
end
a = [];
close all
imshow(frameModel,[]);
hold on
for e = 1:size(pt,1)
    plot(pt(e,2),pt(e,1),'k*');
    text(pt(e,2)+3,pt(e,1),num2str(e));
end
checkBOX = pt([8 6 4 2],:);
qrBOX = pt([3 5 8 7],:);
humanBOX = pt([5 1 7 6],:);
[a(1)] = measurePointBox(checkBOX);
[a(2)] = measurePointBox(qrBOX);
[a(3)] = measurePointBox(humanBOX);
fixedPoints = fliplr(pt);

outSize = size(frameModel);
outP = imref2d(outSize(1:2));

%% collect
CS = [];
for e = 1:numel(qrFileList{1})
    I = imread(qrFileList{1}{e});
    sz = size(I);
    CI = reshape(I,[prod(sz(1:2)) sz(3)]);
    CS = [CS;CI];
    e
end

%%
points = detectMinEigenFeatures(T);
points = detectHarrisFeatures(T,'FilterSize',5);
points = detectSURFFeatures(T);
points = detectKAZEFeatures(T);
[baseFeatures,basePoints] =  extractFeatures(G1,points.Location,'Method','Block','BlockSize',31);
imshow(I,[]);
hold on
plot(basePoints(:,1),basePoints(:,2),'.')
waitforbuttonpress
    
CS = [];
for e = 1:numel(qrFileList{1})
    
    tmpBase = baseFeatures;
    
    I = imread(qrFileList{1}{e});
    G = rgb2gray(I);
    points = detectMinEigenFeatures(G);
    
    sz = size(I);
    CI = reshape(I,[prod(sz(1:2)) sz(3)]);
    T = GMModel.cluster(double(CI));
    T = reshape(T,sz(1:2));
    
    %imshow(T,[]);
    %waitforbuttonpress
    CS = [CS;CI];
    points = detectMinEigenFeatures(T==3);
    
    %{
    Z = squareform(pdist(points.Location));
    pidx = find(sum(Z < 11,1) < 5);
    points = points(pidx);
    %}
    
    
    %[features{e}, valid_points{e},ptVis] = extractHOGFeatures(I,points,'CellSize',[21 21]);
    %[features{e},valid_points{e}] = extractFeatures(G,points);
    
    
    [features{e},valid_points{e}] =  extractFeatures(G,points.Location,'Method','Block','BlockSize',31);
    [indexPairs,matchmetric] = matchFeatures(tmpBase,features{e},'Unique',true,'MaxRatio',.8);
    
    matchedPoints1 = basePoints(indexPairs(:,1),:);
    matchedPoints2 = valid_points{e}(indexPairs(:,2),:);
    %{
    delta = sum((matchedPoints1 - matchedPoints2).^2,2).^.5;
    rm = find(delta > 20);
    
    valid_points{e} = setdiff(valid_points{e},matchedPoints2(rm,:));
    tmpBase = setdiff(tmpBase,matchedPoints1(rm,:));
    
    [indexPairs,matchmetric] = matchFeatures(tmpBase,features{e},'Unique',true,'MaxRatio',.8);
    
    matchedPoints1 = basePoints(kp(indexPairs(:,1)),:);
    matchedPoints2 = valid_points{e}(kp(indexPairs(:,2)),:);
    %}
    close all
    showMatchedFeatures(I1,I,matchedPoints1,matchedPoints2);
    drawnow
    %waitforbuttonpress
    
    
    
    
    %{
    imshow(I,[]);
    hold on
    plot(points.Location(:,1),points.Location(:,2),'.')
    hold off
    drawnow
    %}
end



%%
figure;
flag = true
while flag
    try
        %plot3(CS(1:1000:end,1),CS(1:1000:end,2),CS(1:1000:end,3),'.')
        GMModel = fitgmdist(double(CS(1:100:end,:)),3);
        flag = false;
    catch
        flag = true;
    end
end
%%























