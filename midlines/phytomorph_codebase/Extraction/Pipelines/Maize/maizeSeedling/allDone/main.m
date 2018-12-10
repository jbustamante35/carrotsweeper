clear all
% LINK TO SEARCH FOR - sparkle
%% focus on QR codes
FilePath = '/mnt/tetra/nate/maizeSeedlingData/qrPile/';
FileList = {};
FileExt = {'mat'};
qrPileFileList = gdig(FilePath,FileList,FileExt,1);
for e = 1:numel(qrPileFileList)
    q = load(qrPileFileList{e});
    imshow(q.q.QR,[]);
    drawnow
end
%%
        obj = load(qrPileFileList{4})
        QR = obj.q.QR;
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
%% gather raw NEF from acount
%% gather all sorghumData
dataPath = ['/iplant/home/hirsc213/maizeData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nef')
        wholeFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
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
%FileList = imageFileList;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create arborist for building trees
% init the arborist with rules for building trees
options = statset('Display','iter');
% creaete extract function to give to the arborist
filter = fspecial('gaussian',[21 21],5);
extractFunc = @(X)rgbAndFilterExtract(X,filter);
% create tree growth rules for the arborist
suggestNumberOfClusterFunction = @(data,level,treePara)treePara(1)*(level<=treePara(2)) + 1*(level>treePara(2));
% create cluster parameter generating function for arborist
clusterFunctionFunctionGenerator = @(X,K)fitgmdist(X,K,'Options',options,'Start','plus','RegularizationValue',0.0001);
% feature selection function
idxSelectorFunction = @(X,L)logical([1 1 1 1 1]);
% spec the tree parameters
maxBranch = 3;
maxDepth = 2;
% build the arborist
jA = arborist(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc,maxDepth,maxBranch,idxSelectorFunction);
%% sample data with johnny
SAM = jA.sampleElements(imageFileList,[2 30 4000 .5]);
%% plant forest
forest = jA.plantTrees(SAM,30);
%% tone down filter
filter = fspecial('gaussian',[21 21],1);
extractFunc = @(X)rgbAndFilterExtract(X,filter);
for e = 1:forest.numberOfTrees
    forest.treeSet{e}.extractFunc = extractFunc;
end
%%
close all
for i = 1:30
    oI = double(imread(imageFileList{1}{i}));
    [k] = forest.clusterImage(oI);
    for e = 1:size(k,3)
        rgb = label2rgb(k(:,:,e));
        imshow(cat(2,oI/255,double(rgb)/255),[]);
        drawnow
        waitforbuttonpress
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make whole list
% #3
cnt = 1;
for fl = 1:numel(imageFileList)
    for e = 1:numel(imageFileList{fl})
        wholeList{cnt} = imageFileList{fl}{e};
        cnt = cnt + 1;
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
%% LOOK IN MASS FIX WOW FOR OTHER PARTS NEED TO BREAK DOWN DATASET
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
obj2nd = gmdistribution.fit(CS2(1:100:end,:),6);

%% run GENERATE
parfor R = 1:10000
    simScene(bkFileList,qrFileList,coneFileList,imgSZ,oSIM,R,false);
    R
end
%%
for e = 1:100:numel(cnFileList)
    d = load(cnFileList{e});
    fidx = strfind(cnFileList{e},filesep);
    type = cnFileList{e}(fidx(end-1)+1:fidx(end)-1);
    out = flattenMaskOverlay(d.coneImage,d.coneMask);
    imshow(out,[])
    title(type)
    drawnow
    
    
end
%% load data for cone type network
TEST = false;
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
        if isCone(e)
            tiT = 'without';
        else
            tiI = 'with';
        end
        
        title(tiT);
        drawnow
    end


    e
end
%%
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
%% make histogram for QR for paper
cnt = 1;
HISTO = [];
for e = 1:10:numel(QRFileList)
    try
        d = load(QRFileList{e});
        MSK = imerode(d.MQR,strel('disk',31));
        MSK = bwlarge(MSK);
        R = regionprops(MSK,'BoundingBox');
        mini = imcrop(d.QR,R(1).BoundingBox);
        
        sz = size(mini);
        mmini = reshape(mini,[prod(sz(1:2)) sz(3)]);
        for k = 1:3
            [HISTO(:,cnt,k)] = hist(mmini(:,k),linspace(0,1,256));
        end
        cnt = cnt + 1;
        imshow(mini,[]);
        drawnow
    catch
    end
end
%% make RCN for meta box
% Create image input layer.
inputLayer = imageInputLayer([64 64 3]);
% Define the convolutional layer parameters.
filterSize = [7 7];
numFilters = 10;
% Create the middle layers.
middleLayers = [
    convolution2dLayer(filterSize, numFilters, 'Padding', 1)
    reluLayer()
    convolution2dLayer(filterSize, numFilters, 'Padding', 1)
    reluLayer()
    maxPooling2dLayer(3, 'Stride',2)
    ];

finalLayers = [

    % Add a fully connected layer with 64 output neurons. The output size
    % of this layer will be an array with a length of 64.
    fullyConnectedLayer(64)

    % Add a ReLU non-linearity.
    reluLayer()

    % Add the last fully connected layer. At this point, the network must
    % produce outputs that can be used to measure whether the input image
    % belongs to one of the object classes or background. This measurement
    % is made using the subsequent loss layers.
    fullyConnectedLayer(1+1)

    % Add the softmax loss layer and classification layer.
    softmaxLayer()
    classificationLayer()
];
layers = [
    inputLayer
    middleLayers
    finalLayers
    ];

% Options for step 1.
optionsStage1 = trainingOptions('sgdm', ...
    'MaxEpochs', 10, ...
    'InitialLearnRate', 1e-5, ...
    'CheckpointPath', tempdir,...
    'ExecutionEnvironment','auto','Plots','training-progress');

% Options for step 2.
optionsStage2 = trainingOptions('sgdm', ...
    'MaxEpochs', 10, ...
    'InitialLearnRate', 1e-5, ...
    'CheckpointPath', tempdir,...
    'ExecutionEnvironment','auto','Plots','training-progress');

% Options for step 3.
optionsStage3 = trainingOptions('sgdm', ...
    'MaxEpochs', 10, ...
    'InitialLearnRate', 1e-6, ...
    'CheckpointPath', tempdir,...
    'ExecutionEnvironment','auto','Plots','training-progress');

% Options for step 4.
optionsStage4 = trainingOptions('sgdm', ...
    'MaxEpochs', 10, ...
    'InitialLearnRate', 1e-6, ...
    'CheckpointPath', tempdir,...
    'ExecutionEnvironment','auto','Plots','training-progress');

options = [
    optionsStage1
    optionsStage2
    optionsStage3
    optionsStage4
    ];
%%

%%
detector = trainFasterRCNNObjectDetector(trainTable, layers, options, ...
        'NegativeOverlapRange', [0 0.5], ...
        'PositiveOverlapRange', [.8 1], ...
        'BoxPyramidScale', 4);
    %%
    
options = [
    optionsStage1
    ];
    detector = trainFastRCNNObjectDetector(trainTable, layers, options, ...
        'NegativeOverlapRange', [0 0.5], ...
        'PositiveOverlapRange', [.9 1]);
    %%
       
% Options for step 4.
options = trainingOptions('sgdm', ...
    'MaxEpochs', 10, ...
    'InitialLearnRate', 1e-6, ...
    'CheckpointPath', tempdir,...
    'ExecutionEnvironment','auto');

    detector = trainRCNNObjectDetector(trainTable, layers, options, ...
        'NegativeOverlapRange', [0 0.5], ...
        'PositiveOverlapRange', [.9 1]);
    %% try large
    [FileList] = gdig('/mnt/tetra/nate/AIF/QR/',{},{'tif'},false);
    nT = table;
    for e = 1:numel(FileList)
        I= imread(FileList{e});
        nT(e,'imageFileName') = FileList(e);
        nT(e,'metaBoundingBox') = {{[1 1 size(I,2) size(I,1)]}};
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try large for Faster RCNN - scale the data
close all
[FileList] = gdig('/mnt/tetra/nate/AIF/WHOLE/',{},{'tif'},false);
nT = table;
cnt = 1;
scale = .3;
rScale = '/mnt/tetra/nate/AIF/SCALE/';
mkdir(rScale);
nTF = {};
nTB = {};
for e = 1:numel(FileList)
    %I= imread(FileList{e});
    try
        [~,nm] = fileparts(FileList{e});
        BBfile = strrep(FileList{e},'WHOLE','TAB');
        BBfile = strrep(BBfile,'tif','csv');
        BOX = csvread(BBfile);
     
       
        if scale ~= 1
            I = imread(FileList{e});
            I = imresize(I,scale);
            FileList{e} = [rScale nm '.tif'];
            imwrite(I,FileList{e});
        end
        
        %nT(cnt,'imageFileName') = FileList(e);
        %nT(cnt,'metaDataBox') = {{round([BOX]*scale)}};
       
        nTF{e} = FileList(e);
        nTB{e} = {{round([BOX]*scale)}};
        %{
        I = imread(FileList{cnt});
        dI = insertShape(I, 'Rectangle', round(nT{cnt,'metaDataBox'}{1}),'Color',{'red'});
        imshow(dI,[]);
        drawnow
        cnt = cnt + 1;
        %}
        cnt
        e
    catch ME
        ME
    end
end
%% make table from cells
nT = table;
for e = 1:numel(nTF)
    nT{e,'imageFileName'} = nTF{e};
    nT{e,'metaDataBox'} = {[nTB{1}{1}{1}]};
end
%%
trainTable = nT;
%% train fasterRCNN
trainTable = nT;
detector = trainFasterRCNNObjectDetector(trainTable, layers, options, ...
        'NegativeOverlapRange', [0 0.8], ...
        'PositiveOverlapRange', [.8 1]);
%% try from large format for cascasde
trainCascadeObjectDetector('metaData.xml',nT, ...
'/mnt/tetra/nate/AIF/BIO/','NumCascadeStages',15);
detector = vision.CascadeObjectDetector('metaData.xml');
%% try slow impl and fast train
% Options for step 4.
options = trainingOptions('sgdm', ...
    'MaxEpochs', 10, ...
    'InitialLearnRate', 1e-6, ...
    'CheckpointPath', tempdir,...
    'ExecutionEnvironment','auto');
detector = trainRCNNObjectDetector(nT, layers, options, ...
    'NegativeOverlapRange', [0 0.5], ...
    'PositiveOverlapRange', [.9 1]);
%%
acfDetector = trainACFObjectDetector(nT,'NumStages',20,'NegativeSamplesFactor',100,'ObjectTrainingSize',2*[176 224]);
%% agg - egg
close all;
I = imread(wholeFileList{5200});
%I = imread(trainTable.imageFileName{end});
I = imresize(I,scale);
close all
[BOX,scores] = acfDetector.detect(I);
[scores,sidx] = sort(scores);
detectedImg = insertShape(I, 'Rectangle', BOX,'Color',{'red'});
figure;
imshow(detectedImg);
%%
close all;
I = imread(wholeFileList{300});
%I = imread(trainTable.imageFileName{end});
%I = imresize(I,scale);
I = imresize(I,.3);
%%
close all
[BOX,scores,labels] = detector.detect(I,'MaxSize',[250 250],'NumStrongestRegions',10000);
[scores,sidx] = sort(scores);
detectedImg = insertShape(I, 'Rectangle', BOX,'Color',{'red'});
figure;
imshow(detectedImg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
trainCascadeObjectDetector('stopSignDetector.xml',trainTable(:,[1 2]), ...
'/mnt/tetra/nate/AIF/BIO/','NumCascadeStages',20,'FeatureType','Haar');

detector = vision.CascadeObjectDetector('stopSignDetector.xml');
%%
trainCascadeObjectDetector('stopSignDetector.xml');
%% train table load
tableFile = '/mnt/tetra/nate/maizeSeedlingData/tableFile.mat';
load(tableFile);
%% test crops

I = double(imread(trainTable.wholeImageFileName{1}));
BOX = trainTable.metaDataBox{1};
newSize = 1;
mimcrop(I,BOX,newSize)

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train data for points
close all


fJERK = .85;
JERK = round(40*2.6*fJERK);
disp = false;
%disp = true;
figure;
PStack = {};
near_NStack = {};
far_NStack = {};

masterScale = .5;
scaleJERK = masterScale*JERK;

boxDIM = [31 31];
scale = .1;
sx = round((boxDIM(1)+1)*scale);
sy = round((boxDIM(2)+1)*scale);
grid1 = [];
[grid1(:,:,1),grid1(:,:,2)] = ndgrid(linspace(1,boxDIM(2)+1,sy),linspace(1,boxDIM(1)+1,sx));

parfor e = 1:size(trainTable,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image
    I = double(imread(trainTable.qrImageFileName{e}))/255;
    I = imresize(I,.25);
    
    TALL(:,:,:,e) = imresize(I,.10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the mask for sampling "true" points
    MSK = zeros(size(I,1),size(I,2));
    for pt = 1:size(trainTable.redPoints{e},1)
        MSK(round(.25*trainTable.redPoints{e}(pt,1)-5),round(.25*trainTable.redPoints{e}(pt,2))-5) = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nearPT_RED,farPT_RED,positivePT_RED] = generateSamplePoints_forBOX(MSK,0,.5,10,[10 100],-1);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nearPT_RED] = clearSamplePointsFromBorder(nearPT_RED,size(I),boxDIM,2);
    [farPT_RED] = clearSamplePointsFromBorder(farPT_RED,size(I),boxDIM,2);
    [positivePT_RED] = clearSamplePointsFromBorder(positivePT_RED,size(I),boxDIM,2);
    
    
    
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [nearPT_RED] = point2Box(nearPT_RED,boxDIM);
        [farPT_RED] = point2Box(farPT_RED,boxDIM);
        [positiveBOX_RED] = point2Box(positivePT_RED,boxDIM);
        

        tic
        [PStack{e}] = myCropper(I,positiveBOX_RED,1);
        [near_NStack{e}] = myCropper(I,nearPT_RED,1);
        [far_NStack{e}] = myCropper(I,farPT_RED,1);
        toc
        
        
        if disp
            sI = insertShape(I, 'Rectangle', nearPT_RED,'Color',{'red'},'LineWidth',5);
            sI = insertShape(sI, 'Rectangle', farPT_RED,'Color',{'green'},'LineWidth',5);
            sI = insertShape(sI, 'Rectangle', positiveBOX_RED,'Color',{'blue'},'LineWidth',5);
            imshow(sI,[]);
            drawnow
        end
        e
    catch ME
        ME
        PStack{e} = [];
        near_NStack{e} = [];
        far_NStack{e} = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regress for the QR code center of mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% delta jiggle amount
DELTA = [1000];
% scale range
scaleRange = [.7 1.5];
% over square
OVER_SQUARE = round(2*[1400 1400]);
% resize amount
reSZ = [91 91];
% random samples
RTOT = 10;
[wXi,wYi] = centerSampleJiggle(trainTable,DELTA,scaleRange,OVER_SQUARE,reSZ,RTOT);
%% view for pre-test
close all
for tr = 1:100:size(wXi,4)
    imshow(wXi(:,:,:,tr),[]);
    hold on
    plot(wYi(tr,1)*(OVER_SQUARE(1)/91)^-1+91/2,wYi(tr,2)*(OVER_SQUARE(1)/91)^-1+91/2,'b*')
    drawnow
end
%%
close all
%%%%%%%%%%%%%%%%%%%%%
layers = [imageInputLayer([size(wXi,1) size(wXi,2) size(wXi,3)]);
          convolution2dLayer([11 11],5);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer([11 11],5);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(size(wYi,2));
          regressionLayer();];
options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'InitialLearnRate',.01,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
[nwiY,ZPARA(1,:),ZPARA(2,:)] = zscore(wYi);
qrCenterNet = trainNetwork(wXi,nwiY,layers,options);
%%  test
close all
tr = 800;
for tr = 1:10:size(trainTable,1)
    USE_SQUARE = .9*OVER_SQUARE;
    I = double(imread(trainTable.wholeImageFileName{tr}))/255;
    tic
    sGRID = genIXgrid2(size(I),round(USE_SQUARE/6),[0 0],USE_SQUARE/2);
    
    
    
    [qrLOCS_CORR_P] = fFind_QR(I,qrCenterNet,reSZ,sGRID,USE_SQUARE,ZPARA);
    
    
    
    
    dt = toc;
    tic
    sI = insertShape(I, 'Rectangle', sBOX,'Color',{'red'},'LineWidth',11);
    imshow(sI,[]);
    title(num2str(dt));
    hold on
    plot(qrLOCS_CORR_P(:,1),qrLOCS_CORR_P(:,2),'g*')
    toc
    drawnow
    hold off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regress for the QR code center of mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% delta jiggle amount
DELTA_ZOOM = [500];
% scale range
scaleRange_ZOOM = [.9 1.1];
% over square
OVER_SQUARE_ZOOM = round([1400 1400]);
% resize amount
%reSZ_ZOOM = [201 201];
reSZ_ZOOM = [131 131];
[wXi_ZOOM,wYi_ZOOM] = centerSampleJiggle(trainTable,DELTA_ZOOM,scaleRange_ZOOM,OVER_SQUARE_ZOOM,reSZ_ZOOM,10);

%% regression zoom in
close all
%%%%%%%%%%%%%%%%%%%%%
layers = [imageInputLayer([size(wXi_ZOOM,1) size(wXi_ZOOM,2) size(wXi_ZOOM,3)]);
          convolution2dLayer([11 11],15);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer([11 11],15);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(size(wYi_ZOOM,2));
          regressionLayer();];
options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'InitialLearnRate',.01,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
[NwYi_ZOOM,ZPARA2(1,:),ZPARA2(2,:)] = zscore(wYi_ZOOM);
qrCenterNet2 = trainNetwork(wXi_ZOOM,NwYi_ZOOM,layers,options);
%% cluster
KC = 5;
dY = sum(wYi_ZOOM(:,1:2).^2,2).^.5;
kidx = kmeans(dY,KC);
v = [];
for e = 1:KC
    v(e) = mean(dY(kidx==e));
end
[vs,sidx] = sort(v);
nv = zeros(size(kidx));
for e = 1:KC
    fidx = find(kidx == sidx(e));
    nv(fidx) = e;
end
%% stack extra
%{
sub = 4000;
rnd = randi(size(wXiN,4),[sub 1]);
wXiNT = cat(4,wXi_ZOOM,wXiN(:,:,:,rnd));
nv = [nv ; 6*ones(sub,1)];
%}
%% not used
%[wXiN,wYiN] = NON_centerSampleJiggle(trainTable,DELTA,scaleRange,OVER_SQUARE,reSZ2,totalR);
%%
close all
[wXiN2,wYiN2] = NON_centerSampleJiggle(trainTable,DELTA,scaleRange,OVER_SQUARE,reSZ_ZOOM,5);
%%
sub = 5000;
rnd = randi(size(wXiN2,4),[sub 1]);
wXiNT = cat(4,wXi_ZOOM,wXiN2(:,:,:,rnd));
nvi = [nv ; 6*ones(sub,1)];
%% class zoom in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layers = [imageInputLayer([size(wXiNT,1) size(wXiNT,2) size(wXiNT,3)]);
          convolution2dLayer([11 11],7);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer([11 11],7);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(KC+1);
          softmaxLayer();
          classificationLayer()];
options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',5,...
    'InitialLearnRate',.01,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
minBowlNet = trainNetwork(wXiNT,categorical(nvi),layers,options);
%% view for pre-test
close all
for tr = 1:100:size(wXi,4)
    imshow(wXi(:,:,:,tr),[]);
    hold on
    plot(wYi(tr,1)*(OVER_SQUARE(1)/91)^-1+91/2,wYi(tr,2)*(OVER_SQUARE(1)/91)^-1+91/2,'b*')
    drawnow
end
%% for angle
for e = 1:size(trainTable,1)
    pts = trainTable.redPoints{e};
    vecV1 = pts(3,:) - pts(1,:);
    vecV2 = pts(8,:) - pts(6,:);
    vecH1 = pts(6,:) - pts(1,:);
    vecH2 = pts(8,:) - pts(3,:);
    vecV = mean([vecV1;vecV2]);
    vecH = mean([vecH1;vecH2]);
    vecV(2) = -vecV(2);
    vecH(2) = -vecH(2);
    angleV = atan2(vecV(2),vecV(1))+pi/2;
    angleH = atan2(vecH(2),vecH(1));
    angle = mean([angleV angleH])*180/pi;
    I = imread(trainTable.qrImageFileName{e});
    imshow(I,[]);
    title(num2str(angle));
    drawnow
    A(e) = angle;
end
%%
[nefFileList] = gdig('/mnt/tetra/nate/seedlingDATApile/maizeData/seedlingData/',{},{'nef'},true);

%%  test
close all
tr = 800;


for tr = 1:10:size(trainTable,1)%1:10:numel(nefFileList)%
    
    
    
    USE_SQUARE_ZOOM = OVER_SQUARE_ZOOM;
    USE_SQUARE = .8*OVER_SQUARE;
    %I = double(imread(trainTable.wholeImageFileName{tr}))/255;
    I = double(imread(nefFileList{tr}))/255;
    I(:,(end-100:end),:) = [];
    
    
    VER = randi(size(I,1),1,1);
    HOR = randi(size(I,2),1,1);
    I = circshift(circshift(I,VER,1),HOR,2);
    EXT = USE_SQUARE;
    oI = I;
    I = padarray(I,USE_SQUARE/2,'replicate','both');
    
    tic
    %sGRID = genIXgrid2(size(I),round(USE_SQUARE_ZOOM/6),[0 0],USE_SQUARE_ZOOM/2);
    sGRID = genIXgrid2(size(I),round(USE_SQUARE/4),[0 0],USE_SQUARE/2);
    toc
    
    
    tic
    [qrLOCS_CORR_P_1,scaleOUT_1,sBOX] = fFind_QR(I,qrCenterNet,reSZ,sGRID,USE_SQUARE,ZPARA);
    qrLOCS_CORR_P_1 = flipdim(qrLOCS_CORR_P_1,2);
    fprintf(['Ending first pass scan in:' num2str(toc) '\n']);
    
    
    
    tic
    [qrLOCS_CORR_P,scaleOUT,sBOX,sSAMP] = fFind_QR(I,qrCenterNet2,reSZ_ZOOM,qrLOCS_CORR_P_1,USE_SQUARE_ZOOM,ZPARA2);
    fprintf(['Ending second pass scan in:' num2str(toc) '\n']);
    
    
    tic
    sSAMP = imresize(sSAMP,[201 201]);
    [ptClass,ptProb] = minBowlNet.classify(sSAMP);
    ptValue = ptProb*(1:6)';
    [sV,sidx] = sort(ptValue);
    sidx = sidx(3);
    fprintf(['Ending classification pass scan in:' num2str(toc) '\n']);
    
    
    
    
    
    
    qrLOCS_CORR_P = bsxfun(@minus,qrLOCS_CORR_P,USE_SQUARE/2);
    qrLOCS_CORR_P_1 = bsxfun(@minus,qrLOCS_CORR_P_1,USE_SQUARE/2);
    
    
    selectBOX = point2Box(qrLOCS_CORR_P(sidx,:),USE_SQUARE_ZOOM*scaleOUT(sidx));
    
    
    redBOX = point2Box(qrLOCS_CORR_P(sidx,:),boxDIM);
    qrForRed = mimcrop(oI,redBOX,.2);
    redPT = rcNet.predict(qrForRed);
    redPT = bsxfun(@times,redPT,redPARA2(2,:));
    redPT = bsxfun(@plus,redPT,redPARA2(1,:));
    
    
    imshow(qrForRed,[]);
    hold on;
    %sredPT = bsxfun(@minus,sredPT,redBOX(1:2));
    sredPT = redPT*.2;
    sredPT = reshape(sredPT,[8 2]);
    plot(sredPT(:,1),sredPT(:,2),'r.')
    redPT = reshape(redPT,[8 2]);
    redPT = bsxfun(@plus,redPT,selectBOX(1:2));
    
    
    imshow(oI,[]);
    hold on
    plot(redPT(:,1),redPT(:,2),'r.');
    hold off
    
    
    tic
    sI = insertShape(oI, 'Rectangle', selectBOX,'Color',{'red'},'LineWidth',11);
    imshow(sI,[]);
    title(num2str(dt));
    hold on
    plot(qrLOCS_CORR_P(sidx,1),qrLOCS_CORR_P(sidx,2),'g*')
    plot(qrLOCS_CORR_P_1(sidx,2),qrLOCS_CORR_P_1(sidx,1),'r*')
    
    
    plot(sGRID(:,2),sGRID(:,1),'k.');
    plot(qrLOCS_CORR_P_1(:,2),qrLOCS_CORR_P_1(:,1),'k.');
    plot(qrLOCS_CORR_P(:,1),qrLOCS_CORR_P(:,2),'k.')
    
   
    for e = 1:size(sGRID,1)
        plot([(sGRID(e,2)) qrLOCS_CORR_P_1(e,2)],[(sGRID(e,1)) qrLOCS_CORR_P_1(e,1)],'k')
        plot([(qrLOCS_CORR_P_1(e,2)) qrLOCS_CORR_P(e,1)],[(qrLOCS_CORR_P_1(e,1)) qrLOCS_CORR_P(e,2)],'k')
    end
    
    
    
    e = sidx;
    plot([(sGRID(e,2)) qrLOCS_CORR_P_1(e,2)],[(sGRID(e,1)) qrLOCS_CORR_P_1(e,1)],'g')
    plot([(qrLOCS_CORR_P_1(e,2)) qrLOCS_CORR_P(e,1)],[(qrLOCS_CORR_P_1(e,1)) qrLOCS_CORR_P(e,2)],'g')
    
    toc
    drawnow
    hold off
    %waitforbuttonpress
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% red corner jiggle - find the corners - red poimt regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
jiggelVec = [-150 150];
boxDIM = [1130 937];
parfor e = 1:size(trainTable,1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image
    I = double(imread(trainTable.wholeImageFileName{e}))/255;
    RP = trainTable.redPoints{e};
    RTMP_X = [];
    RTMP_Y = [];
    for r = 1:10
        JG = jiggelVec(2)*(rand([1 2]) - .5);
        rBOX = trainTable.metaDataBox{e};
        rBOX(3:4) = (boxDIM);
        rBOX(1:2) = rBOX(1:2) + JG;
        
        % red points
        tRP = bsxfun(@minus,RP,round(JG))-5;
        %tRP = RP;
        
        pts = trainTable.redPoints{e};
        vecV1 = pts(3,:) - pts(1,:);
        vecV2 = pts(8,:) - pts(6,:);
        vecH1 = pts(6,:) - pts(1,:);
        vecH2 = pts(8,:) - pts(3,:);
        vecV = mean([vecV1;vecV2]);
        vecH = mean([vecH1;vecH2]);
        vecV(2) = -vecV(2);
        vecH(2) = -vecH(2);
        angleV = atan2(vecV(2),vecV(1))+pi/2;
        angleH = atan2(vecH(2),vecH(1));
        angle = mean([angleV angleH])*180/pi;
        CP = mean(pts,1);
        tRP = bsxfun(@minus,tRP,CP);
        tRP = bsxfun(@plus,tRP,boxDIM/2);
        
        
        CP = CP + round(JG)-5;
        CP = CP + trainTable.metaDataBox{e}(1:2);
        
        redBOX = point2Box(CP,boxDIM);
        
       
        %tRP = bsxfun(@plus,tRP,trainTable.metaDataBox{e}(1:2));
        
        
        tmp = mimcrop(I,redBOX,.2);
        %imshow(tmp,[]);
        rBOX;
        
        
        %{
        imshow(tmp,[]);
        hold on
        plot(tRP(:,1),tRP(:,2),'g*')
        drawnow
        hold off
        %}
        size(tmp);
        %tmp = imresize(tmp,.1);
        
        
        RTMP_X(:,:,:,r) = tmp;
        RTMP_Y(:,:,r) = tRP;
    end
    XX{e} = RTMP_X;
    YY{e} = RTMP_Y;
    e
end
%
sz = size(XX{1});
per = sz(4);
sz(4) = sz(4)*numel(XX);
rX = zeros(sz);
rY = [];
for e = 1:numel(XX)
    str = (e-1)*per + 1;
    stp = str + per - 1;
    rX(:,:,:,str:stp) = XX{e};
    tmp = YY{e};
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);
    rY = [rY;tmp'];
    e
end
%{
%%
close all
tr = 200;
for tr = 1:10:size(rX,4)
    imshow(rX(:,:,:,tr),[]);
    tmpPR =.2* reshape(rY(tr,:),[8 2]);
    hold on
    plot(tmpPR(:,1),tmpPR(:,2),'r.')
    drawnow
    hold off
end
%}
%%
for e = 1:size(rX,4)
    srX(:,:,:,e) = imresize(rX(:,:,:,e),.5);
    for k = 1:3
        tmp = srX(:,:,k,e);
        srX(:,:,k,e) = srX(:,:,k,e) - mean(tmp(:));
    end
end
%%
SKIP = 10;
TR = size(rX,4)/SKIP;
wSZ = [13 13];
nSZ = [(size(rX,1)-(wSZ(1)-1)) (size(rX,2)-(wSZ(2)-1))];
D = zeros(TR,prod(nSZ),prod(wSZ));
STACK = zeros(size(P,2),prod(wSZ)^2,TR);
cnt = 1;
for tr = 1:SKIP:TR*SKIP
    tmp = im2col(rX(:,:,1,tr),wSZ,'sliding')';
    R = zeros(size(tmp,1),size(tmp,2)^2);
    tic
    for p = 1:size(tmp,1)
        tempy = tmp(p,:)'*tmp(p,:);
        R(p,:) = tempy(:)';
    end
    toc
    STACK(:,:,cnt)  = mtimesx(P,'T',R);
    cnt = cnt + 1;
    tr
end
%% logic gate
wSZ = [13 13];
SKIP = 10;
TR = size(rX,4)/SKIP;
nSZ = [(size(rX,1)-(wSZ(1)-1)) (size(rX,2)-(wSZ(2)-1))];
logicData = zeros(prod(wSZ),TR,'uint8');
targetLogic = zeros(TR,1,'uint8');
str = 1;
BLOCK = prod(nSZ);
cnt = 1;
for tr = 1:SKIP:TR*SKIP
    
    tmpPR = .2*reshape(rY(tr,:),[8 2]);
    
    %{
    imshow(rX(:,:,1,tr),[]);
    hold on
    plot(tmpPR(1,1),tmpPR(1,2),'r*');
    drawnow
    %}
    
    
    msk = zeros(size(rX,1),size(rX,2));
    msk(round(tmpPR(1,2)),round(tmpPR(1,1))) = 1;
    msk = im2col(msk,wSZ,'sliding');
    msk = imdilate(msk,strel('disk',2,0));
    msk = msk((end-1)/2,:);
    
    
    stp = str + BLOCK - 1;
    tmp = im2col(255*rX(:,:,1,tr),wSZ,'sliding');
    logicData(:,str:stp) = tmp;
    targetLogic(str:stp) = msk(:);
    str = stp + 1;
    cnt = cnt + 1
end
logicData = logicData';
%% generate random gate
gate = randi(255,[1 size(logicData,2)],'uint8');
%%
tic
[res] = logicFunc(gate,logicData,targetLogic);
toc
%% stack lew logical
lX = zeros(size(logicData,1),size(logicData,2)*8,'logical');
cnt = 1;
for pos = 1:size(logicData,2)
    for m = 0:7
        lX(:,cnt) = bitget(logicData(:,pos),m+1);
        cnt = cnt + 1
    end
end
%%
cnt = 1;
lXtmp = zeros(size(tmp,1),size(tmp,2)*8,'logical');
for pos = 1:size(tmp,2)
    for m = 0:7
        lXtmp(:,cnt) = bitget(uint8(tmp(:,pos)),m+1);
        cnt = cnt + 1
    end
 end
%%
IDX1 = find(targetLogic);
IDX0 = find(targetLogic==0);
MAG = 5;
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX1;IDX0(1:(MAG*numel(IDX1)))];
%%
%[IDXS, ZS] = rankfeatures(double(logicData(IDX,:))',logical(targetLogic(IDX))','Criterion','entropy');
%[IDXS, ZS] = rankfeatures(lX',logical(targetLogic(IDX))','Criterion','entropy');
for e = 1:size(lX,2)
    ME(e) = mutualinfo(uint8(targetLogic),uint8(lX(:,e)));
    e
end
%%
%%[IDX, Z] = rankfeatures(lX, targetLogic,'Criterion','entropy');
%%
%mdl = fitglm(double(logicData(IDX,IDXS(1:2))),double(targetLogic(IDX)),'Link','probit','Distribution','binomial');
%mdl = stepwiseglm(double(logicData(IDX,IDXS(1:10))),logical(targetLogic(IDX)),'y~1','Distribution','binomial','upper','poly2222222222');
%%
tmpX = lX(IDX,find(inmodel));
mdl = stepwiseglm(tmpX,logical(targetLogic(IDX)),'y~1','Distribution','binomial','upper','quadratic');
%%
OP = statset('UseParallel',true,'Display','iter');
[inmodel,history] = sequentialfs(@mySelectionProcess,lX(IDX,:),targetLogic(IDX),'options',OP,'cv',10,'NFeatures',10);
%%
selV = zeros(1,size(lX,2),'logical');
for r = 1:10
    EN = [];
    
    for e = 1:size(lX,2)
        selVT = selV;
        selVT(e) = true;
        EN(e) = mySelectionProcess(lX(IDX,selVT),targetLogic(IDX));
        plot(EN)
        drawnow
    end
    
    [~,midx] = min(EN);
    selV(midx) = true;
end

%%
func = @(g)logicFunc(g,logicData,targetLogic);
options = optimoptions('particleswarm','UseParallel',true,'Display','iter','SwarmSize',100);
[x,fval,exitflag] = particleswarm(func,169,zeros(1,size(logicData,2)),255*ones(1,size(logicData,2)),options);
%%
A = squeeze(STACK(1,:,:))';
B = rY(1:10:end,9:end);
%%
options = statset('UseParallel',true,'Display','Iter');
%%
%%
%%
[nY,nnU,nnS]= zscore(rY);
%%
B = nY(1:10:end,1:8);
Ax = squeeze(STACK(2,:,:))';
[XL,YL,XS,YS,BETAx] = plsregress(Ax,B,15);
[Lx,FitInfoX] = lasso(Ax,B(:,1),'Options',options,'Alpha',.0001,'CV',10);
[~,Cx,Ux,Ex] = PCA_FIT_FULL(Ax,50);
regRx = Cx\B;



B = nY(1:10:end,9:end);
Ay = squeeze(STACK(1,:,:))';
[XL,YL,XS,YS,BETAy] = plsregress(Ay,B,15);
[Ly,FitInfoY] = lasso(Ay,B,'Options',options,'Alpha',.0001,'CV',10);
[~,Cy,Uy,Ey] = PCA_FIT_FULL(Ay,50);
regRy = Cy\B;
%%
for tr = 1:10:size(rX,4)
    toM = rX(:,:,1,tr);
    reY = rY(tr,:);
    tmp = im2col(toM,wSZ,'sliding')';
    
    %tmp = im2col(255*rX(:,:,1,tr),wSZ,'sliding');
    
    R = zeros(size(tmp,1),size(tmp,2)^2);
    tic
    for p = 1:size(tmp,1)
        tempy = tmp(p,:)'*tmp(p,:);
        R(p,:) = tempy(:)';
    end

    %{
    R = PCA_REPROJ(R,E,U);
    rho = mtimesx(R,xR);


    RyC = PCA_REPROJ(R,Ey,Uy);
    yI = RyC*regRy;
    for e = 1:size(yI,2)
        Iyy(:,:,e) = col2im(yI(:,e),wSZ,size(toM),'sliding');
    end
    RxC = PCA_REPROJ(R(2,:),Ex,Ux);
    %}

    %
    R = mtimesx(P,'T',R);
    %
    iX = [1 R(2,:)]*BETAx;
    iY = [1 R(1,:)]*BETAy;
    iPT = [iX iY];
    Ry = PCA_REPROJ(R(1,:),Ey,Uy);
    Rx = PCA_REPROJ(R(2,:),Ex,Ux);


    PT = [(Rx*regRx)' (Ry*regRy)'];
    PT = PT(:)';
    PT = bsxfun(@times,PT,nnS);
    PT = bsxfun(@plus,PT,nnU);

    iPT = bsxfun(@times,iPT,nnS);
    iPT = bsxfun(@plus,iPT,nnU);


    iPT = reshape(iPT',[8 2]);
    PT = reshape(PT',[8 2]);
    PTU = reshape(nnU',[8 2]);
    reY = reshape(reY',[8 2]);
    
    PT = PT * .2*2;
    iPT = iPT * .2*2;
    reY = reY * .2*2;
    PTU = PTU * .2*2;
    SH = imresize(rX(:,:,:,tr),2);
    close all
    imshow(SH,[]);
    hold on
    plot(PT(:,1),PT(:,2),'g.');
    plot(iPT(:,1),iPT(:,2),'k.');
    plot(reY(:,1),reY(:,2),'b.');
    %plot(PTU(:,1),PTU(:,2),'r.');
    drawnow
end
%%
pSTORE = PT;
%%

plot(pSTORE(:,1),pSTORE(:,2),'ro');
%%
x = optimvar(x,13,1)
%%
[iG(:,:,1),iG(:,:,2)] = ndgrid((wSZ(1)-1)/2+1:(size(rX,1)-(wSZ(1)-1)/2),(wSZ(2)-1)/2+1:(size(rX,2)-(wSZ(2)-1)/2));
P = [reshape(iG(:,:,1),[size(iG,1)*size(iG,2) 1]) reshape(iG(:,:,2),[size(iG,1)*size(iG,2) 1])];
tmp = mtimesx(X,V);


oldSchool(D,ones(13*13,5),P);
%%
%%%%%%%%%%%%%%%%%%%%%
layers = [imageInputLayer([size(srX,1) size(srX,2) size(srX,3)],'Normalization','none');
          convolution2dLayer([11 11],5);
          reluLayer();
          fullyConnectedLayer(size(rY,2));
          regressionLayer();];
options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.15,...
    'InitialLearnRate',.0000001,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
[redY,redPARA2(1,:),redPARA2(2,:)] = zscore(rY);
rcNet = trainNetwork(srX,redY,layers,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stack and train for red corner network
G1D1 = [];
for e = 1:numel(PStack)
    if ~isempty(PStack{e})
        if size(PStack{e}{1},4) ~= 8
            e
        end
        G1D1 = cat(4,G1D1,PStack{e}{1});
    end
end
G2D1 = [];
for e = 1:numel(near_NStack)
    if ~isempty(near_NStack{e})
        G2D1 = cat(4,G2D1,near_NStack{e}{1});
    end
end
G3D1 = [];
for e = 1:numel(far_NStack)
    if ~isempty(far_NStack{e})
        G3D1 = cat(4,G3D1,far_NStack{e}{1});
    end
end

G1D2 = [];
for e = 1:numel(PStack)
    if ~isempty(PStack{e})
        G1D2 = cat(4,G1D2,PStack{e}{2});
    end
end
G2D2 = [];
for e = 1:numel(near_NStack)
    if ~isempty(near_NStack{e})
        G2D2 = cat(4,G2D2,near_NStack{e}{2});
    end
end
G3D2 = [];
for e = 1:numel(far_NStack)
    if ~isempty(far_NStack{e})
        G3D2 = cat(4,G3D2,far_NStack{e}{2});
    end
end
%{
G1D3 = [];
for e = 1:numel(PStack)
    if ~isempty(PStack{e})
        G1D3 = cat(4,G1D3,PStack{e}{3});
    end
end
G2D3 = [];
for e = 1:numel(near_NStack)
    if ~isempty(near_NStack{e})
        G2D3 = cat(4,G2D3,near_NStack{e}{3});
    end
end
G3D3 = [];
for e = 1:numel(far_NStack)
    if ~isempty(far_NStack{e})
        G3D3 = cat(4,G3D3,far_NStack{e}{3});
    end
end
%}
G3D3 = [];
%Y = [2*ones(size(G1D1,4),1);1*ones(size(G2D1,4),1);0*ones(size(G3D1,4),1)];
Y = [repmat((1:8)',[size(G1D1,4)/8 1]);9*ones(size(G2D1,4),1);10*ones(size(G3D1,4),1)];
%X = cat(3,cat(4,G1D1,G2D1,G3D1),cat(4,G1D2,G2D2,G3D2),cat(4,G1D3,G2D3,G3D3));
X = cat(3,cat(4,G1D1,G2D1,G3D1),cat(4,G1D2,G2D2,G3D2));

layers = [imageInputLayer([size(X,1) size(X,2) size(X,3)]);
          convolution2dLayer([7 7],15);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(10);
          softmaxLayer();
          classificationLayer()];
options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'InitialLearnRate',.0001,...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
cNet = trainNetwork(X,categorical(Y),layers,options);
%% test network
close all
boxDIM = [31 31];
init_QR_spacing = 10;
I = double(imread(trainTable.qrImageFileName{400}))/255;
I = imresize(I,.25);
[g1,g2] = ndgrid(1:size(I,1),1:size(I,2));
g = mod(g1,init_QR_spacing) == 0 & mod(g2,init_QR_spacing) == 0;

g(1:(ceil(boxDIM(2)/2)+1),:) = 0;
g((end-ceil(boxDIM(2)/2)-1):end,:) = 0;

g(:,1:(ceil(boxDIM(1)/2)+1)) = 0;
g(:,(end-ceil(boxDIM(1)/2)-1):end) = 0;

gX = [];
[gX(:,2),gX(:,1)] = find(g);
newSZ = [sum(any(g,2)) sum(any(g,1))];
[gBOX] = point2Box(gX,boxDIM);
featureFunction = @(X,n)featureFunction_wt(X,n,1,3,'haar',false);
tic
[gSTACK] = myCropper(I,gBOX,1,true);
toc
%%
%MS = cat(3,gSTACK{1},gSTACK{2},gSTACK{3});
MS = cat(3,gSTACK{1},gSTACK{2});
[ptL,ptP] = cNet.classify(MS);
K = reshape(ptL,newSZ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather data for scale of object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScaleMax = 2;
ScaleMin = .3;
ScaleRange = ScaleMax - ScaleMin;
ScaleMean = .5*(ScaleMax + ScaleMin);

boxDIM = [1030 837];
X = {};
Y = {};
parfor e = 1:size(trainTable,1)
    
    I = double(imread(trainTable.wholeImageFileName{e}))/255;
    sampleBox = trainTable.metaDataBox{e};
    
    
    oS = [];
    dis = [];
    OFFSET_DIAG = [];
    for r = 1:10
        rScale = (rand(1)-.5)*ScaleRange + ScaleMean;
        centerPoint = round(mean([sampleBox(1:2);sampleBox(1:2)+sampleBox(3:4)]));
        tmpDIM = boxDIM*rScale^-1;


        cross1 = norm(trainTable.redPoints{e}(1,:) - trainTable.redPoints{e}(8,:));
        cross2 = norm(trainTable.redPoints{e}(3,:) - trainTable.redPoints{e}(6,:));
        dis(r) = .5*(cross1 + cross2)*rScale;
        
        

        
        

        [scaleBOX] = point2Box(centerPoint,tmpDIM);

        OFFSET_DIAG(r) = norm(scaleBOX(3:4)) - dis(r);


        tmpI = imcrop(I,scaleBOX);
        tmpI = imresize(tmpI,flipdim(boxDIM,2));
        oS(:,:,:,r) = imresize(tmpI,.1);
    end
    %imshow(tmpI,[]);
    %drawnow
    
    Y{e} = dis;
    YP{e} = OFFSET_DIAG;
    X{e} = oS;
    
    e
end
%%
iX = [];
iY = [];
ipY = [];
for e = 1:numel(X)
   iX =  cat(4,iX,X{e});
   iY = [iY,Y{e}];
   ipY = [ipY,YP{e}];
   e
end
mean(ipY)
%%
%%%%%%%%%%%%%%%%%%%%%
layers = [imageInputLayer([size(iX,1) size(iX,2) size(iX,3)]);
          convolution2dLayer([11 11],7);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer([11 11],3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(1);
          regressionLayer();];
options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'InitialLearnRate',.000001,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
scaleNet = trainNetwork(iX,(iY-mean(iY))',layers,options);
%% sim scale net
tic
tr = 900;
scaleValue = double(scaleNet.predict(iX(:,:,:,tr)))+mean(iY);
toc
delta = scaleValue - iY(tr);
percentDelta = abs(delta)*iY(tr)^-1;
delta
percentDelta
%% test all
scaleValue = double(scaleNet.predict(iX));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train data for QR detection
close all
fJERK = .85;
JERK = round(40*2.6*fJERK);
disp = false;
featureFunction = @(X,n)featureFunction_wt(X,n,.1,3,'haar',disp);
figure;
PStack = {};
near_NStack = {};
far_NStack = {};

masterScale = .5;
scaleJERK = masterScale*JERK;

boxDIM = [1030 837];
scale = .1;
sx = round((boxDIM(1)+1)*scale);
sy = round((boxDIM(2)+1)*scale);
grid1 = [];
[grid1(:,:,1),grid1(:,:,2)] = ndgrid(linspace(1,boxDIM(2)+1,sy),linspace(1,boxDIM(1)+1,sx));

parfor e = 1:size(trainTable,1)
    
    I = double(imread(trainTable.wholeImageFileName{e}))/255;
    %I = imresize(I,masterScale);
    
    sampleBox = trainTable.metaDataBox{e};
    boxDIM = sampleBox(3:4);
    boxDIM = [1030 837];
    centerPoint = round(mean([sampleBox(1:2);sampleBox(1:2)+sampleBox(3:4)]));
    
    pointMask = zeros(size(I,1),size(I,2));
    pointMask(centerPoint(2),centerPoint(1)) = 1;
    [nearPT,farPT,positivePT] = generateSamplePoints_forBOX(pointMask,JERK,.25,300,[10 10],40);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rm = nearPT(:,2) >= (size(I,1) - boxDIM(2)/2-2);
    nearPT(rm,:) = [];
    
    rm = nearPT(:,1) >= (size(I,2) - boxDIM(1)/2-2);
    nearPT(rm,:) = [];
    
    rm = nearPT(:,2) <= boxDIM(2)/2+2;
    nearPT(rm,:) = [];
    
    rm = nearPT(:,1) <= boxDIM(1)/2+2;
    nearPT(rm,:) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rm = farPT(:,2) >= (size(I,1) - boxDIM(2)/2-2);
    farPT(rm,:) = [];
    
    rm = farPT(:,1) >= (size(I,2) - boxDIM(1)/2-2);
    farPT(rm,:) = [];
    
    rm = farPT(:,2) <= boxDIM(2)/2+2;
    farPT(rm,:) = [];
    
    rm = farPT(:,1) <= boxDIM(1)/2+2;
    farPT(rm,:) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rm = positivePT(:,2) >= (size(I,1) - boxDIM(2)/2-2);
    positivePT(rm,:) = [];
    
    rm = positivePT(:,1) >= (size(I,2) - boxDIM(1)/2-2);
    positivePT(rm,:) = [];
    
    rm = positivePT(:,2) <= boxDIM(2)/2+2;
    positivePT(rm,:) = [];
    
    rm = positivePT(:,1) <= boxDIM(1)/2+2;
    positivePT(rm,:) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    try
        [nearPT] = point2Box(nearPT,boxDIM);
        [farPT] = point2Box(farPT,boxDIM);
        [positiveBOX] = point2Box(positivePT,boxDIM);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tic
        [PStack{e}] = myCropper(I,positiveBOX,.1);
        [near_NStack{e}] = myCropper(I,nearPT,.1);
        [far_NStack{e}] = myCropper(I,farPT,.1);
        toc
        
        
        if disp
            sI = insertShape(I, 'Rectangle', falseBOX,'Color',{'red'},'LineWidth',11);
            sI = insertShape(sI, 'Rectangle', positiveBOX,'Color',{'green'},'LineWidth',11);
            imshow(sI,[]);
            drawnow
        end
        e
    catch ME
        ME
        PStack{e} = [];
        near_NStack{e} = [];
        far_NStack{e} = [];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stack and train for QR network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1D1 = [];
for e = 1:numel(PStack)
    if ~isempty(PStack{e})
        G1D1 = cat(4,G1D1,PStack{e}{1});
    end
end
G2D1 = [];
for e = 1:numel(near_NStack)
    if ~isempty(near_NStack{e})
        G2D1 = cat(4,G2D1,near_NStack{e}{1});
    end
end
G3D1 = [];
for e = 1:numel(far_NStack)
    if ~isempty(far_NStack{e})
        G3D1 = cat(4,G3D1,far_NStack{e}{1});
    end
end

G1D2 = [];
for e = 1:numel(PStack)
    if ~isempty(PStack{e})
        G1D2 = cat(4,G1D2,PStack{e}{2});
    end
end
G2D2 = [];
for e = 1:numel(near_NStack)
    if ~isempty(near_NStack{e})
        G2D2 = cat(4,G2D2,near_NStack{e}{2});
    end
end
G3D2 = [];
for e = 1:numel(far_NStack)
    if ~isempty(far_NStack{e})
        G3D2 = cat(4,G3D2,far_NStack{e}{2});
    end
end
%{
G1D3 = [];
for e = 1:numel(PStack)
    if ~isempty(PStack{e})
        G1D3 = cat(4,G1D3,PStack{e}{3});
    end
end
G2D3 = [];
for e = 1:numel(near_NStack)
    if ~isempty(near_NStack{e})
        G2D3 = cat(4,G2D3,near_NStack{e}{3});
    end
end
G3D3 = [];
for e = 1:numel(far_NStack)
    if ~isempty(far_NStack{e})
        G3D3 = cat(4,G3D3,far_NStack{e}{3});
    end
end
%}
G3D3 = [];
Y = [2*ones(size(G1D1,4),1);1*ones(size(G2D1,4),1);0*ones(size(G3D1,4),1)];
%X = cat(3,cat(4,G1D1,G2D1,G3D1),cat(4,G1D2,G2D2,G3D2),cat(4,G1D3,G2D3,G3D3));
X = cat(3,cat(4,G1D1,G2D1,G3D1),cat(4,G1D2,G2D2,G3D2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layers = [imageInputLayer([size(X,1) size(X,2) size(X,3)]);
          convolution2dLayer([7 7],15);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(3);
          softmaxLayer();
          classificationLayer()];
options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
qrNet = trainNetwork(X,categorical(Y),layers,options);
%% test network single
close all
I = double(imread(trainTable.wholeImageFileName{400}))/255;
[g1,g2] = ndgrid(1:size(I,1),1:size(I,2));
g = mod(g1,2*JERK) == 0 & mod(g2,2*JERK) == 0;

g(1:(ceil(boxDIM(2)/2)+1),:) = 0;
g((end-ceil(boxDIM(2)/2)-1):end,:) = 0;

g(:,1:(ceil(boxDIM(1)/2)+1)) = 0;
g(:,(end-ceil(boxDIM(1)/2)-1):end) = 0;

gX = [];
[gX(:,2),gX(:,1)] = find(g);
newSZ = [sum(any(g,2)) sum(any(g,1))];
[gBOX] = point2Box(gX,boxDIM);
featureFunction = @(X,n)featureFunction_wt(X,n,.1,3,'haar',false);
tic
[gSTACK] = myCropper(I,gBOX,.1);
toc
%MS = cat(3,gSTACK{1},gSTACK{2},gSTACK{3});
MS = cat(3,gSTACK{1},gSTACK{2});
%% test netowrk all
[pY,ppY] = qrNet.classify(X);
pppY = double(ppY)*(1:3)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% best test
close all
I = double(imread(trainTable.wholeImageFileName{30}))/255;
I = imresize(I,.5);
[dx,dy,V] = impixel(I);
%%
close all
pt = [1 dx dy];
MAG = [2 4000 4000];
func = @(X)myQRgrade(X,MAG,I,[1030 837],pt,{qrNet,scaleNet,rcNet},mean(iY),[mean(iY) std(iY)]);
options = optimset('Display','iter','MaxIter',100);
[xo,~,~,~,STORE] = fminsearch(func,MAG,options);

[~,parts,value,RP] = func(xo);


RP = RP + reshape(urY,[8 2]);


xoT = bsxfun(@plus,xo,pt - MAG);
xoBOX = [xoT(:,2:3) [[1030 837]'*(xoT(:,1))']'];
RP = RP * (xoT(:,1));
xoBOX = point2Box(xoBOX(1:2),xoBOX(3:4));
dI = insertShape(I, 'Rectangle', xoBOX,'Color',{'red'},'LineWidth',11);

RP = bsxfun(@plus,RP,xoBOX(:,1:2));
imshow(dI,[]);
hold on
plot(RP(:,1),RP(:,2),'b*');
plot(pt(2),pt(3),'g*')
%% PSO
close all
pt = [1 dx dy];
MAG = [6 4000 4000];
func = @(X)myQRgrade(X,MAG,I,[1030 837],pt,{qrNet,scaleNet,rcNet},mean(iY),[mean(iY) std(iY)]);


options = optimoptions('particleswarm','SwarmSize',30,'Display','iter','MaxIterations',35,'UseParallel',false,'UseVectorized',true);
xo = particleswarm(func,3,[5.5;3500;3500],[8;4500;4500],options);
[~,parts,value,RP] = func(xo);

RP = RP + reshape(urY,[8 2]);



xoT = bsxfun(@plus,xo,pt - MAG);
xoBOX = [xoT(:,2:3) [[1030 837]'*(xoT(:,1))']'];
xoBOX = point2Box(xoBOX(1:2),xoBOX(3:4));
dI = insertShape(I, 'Rectangle', xoBOX,'Color',{'red'},'LineWidth',11);
imshow(dI,[]);
hold on
RP = RP * (xoT(:,1));
RP = bsxfun(@plus,RP,xoBOX(:,1:2));
plot(pt(2),pt(3),'g*')
plot(RP(:,1),RP(:,2),'b*');
%%
close all
for e = 1:1:size(STORE,3)
    tmpS = STORE(:,:,e)';
    tmpS = bsxfun(@plus,tmpS,pt - MAG);
    
    
    xoBOX = [xoT(:,2:3) [[1030 837]'*(xoT(:,1))']'];
    
    for b = 1:size(dBOX,1)
        dBOX(b,:) = point2Box(dBOX(b,1:2),dBOX(b,3:4));
    end
    
    
    dI = insertShape(I, 'Rectangle', dBOX,'Color',{'red'},'LineWidth',11);
    imshow(dI,[]);
    drawnow
    pause(.1)
end


%%
[detector] = myRCNtrain(trainTable,STACK,'./output/');
%%
func = cFlow('myRCNtrain');
func(trainTable,STACK,'./output/');
func.isGPU(1);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(100,100);
    %%
    close all
    for e = 1:size(trainTable,1)
        I = imread(trainTable.imageFileName{e});
        STACK(:,:,:,e) = I;
        dI = insertShape(I, 'Rectangle', trainTable.metaDataBox{e},'Color',{'red'});
        dI = insertShape(dI, 'Rectangle', trainTable.staticDataBox{e},'Color',{'green'});
        dI = insertShape(dI, 'Rectangle', trainTable.dynamicDataBox{e},'Color',{'blue'});
        imshow(dI,[]);
        drawnow
    end
    %% 
    close all
    
    %I = imread(trainTable.imageFileName{1});
    I = imread(wholeFileList{10000});
    %I = imread(trainTable.imageFileName{2});
    %I = imread(FileList{e});
    %%
    rI = imread(nT.imageFileName{1});
    %rI = I;
    close all
    %rI = imresize(rI,10);
    %detector.MergeThreshold = 4;
    [BOX] = detector.step(rI);
    detectedImg = insertShape(rI, 'Rectangle', BOX,'Color',{'r'});
    figure
    imshow(detectedImg)
    %%
    close all;
    I = imread(wholeFileList{300});
    %I = imread(trainTable.imageFileName{end});
    I = imresize(I,scale);
    [BOX,scores,labels] = detector.detect(I);
    [scores,sidx] = sort(scores);
    detectedImg = insertShape(I, 'Rectangle', BOX(sidx(end),:),'Color',{'red'});
    figure;
    imshow(detectedImg);
%%
IDX = randi([1 numel(wholeFileList)],[1 5000]);
rndList = wholeFileList(IDX);
%%
I = imread(nefFileList{1});
%%
cwtstruct = cwtft2(double(I)/255);
%%
dwtmode('sp0');
[C,S] = wavedec2(double(I)/255,3,'haar');
%%
[H1,V1,D1] = detcoef2('all',C,S,3);
A1 = appcoef2(C,S,'haar',2);
%% play
I = double(imread(nefFileList{1}))/255;
%% local focus
FilePath = '/mnt/tetra/nate/seedlingDATApile/maizeData/seedlingData/';
FileList = {};
FileExt = {'nef'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
dwtmode('sp0');
tic
L = 1;
[C,S] = wavedec2(double(I)/255,L,'haar');
[H1,V1,D1] = detcoef2('all',C,S,L);
toc

%%
[x,y,v] = impixel(I);
pointMask = zeros(size(I,1),size(I,2));
pointMask(y,x) = 1;
[falsePT,positivePT] = generateSamplePoints_forBOX(pointMask,21,.25,300,[10 10],10);
[falseBOX] = point2Box(falsePT,[100 100]);
[positiveBOX] = point2Box(positivePT,[100 100]);
sI = insertShape(I, 'Rectangle', falseBOX,'Color',{'red'},'LineWidth',3);
sI = insertShape(sI, 'Rectangle', positiveBOX,'Color',{'green'},'LineWidth',3);
imshow(sI,[]);
%%
featureFunction = @(I)featureFunction_wt(I,.25,1,'haar');
BOX = [300 400 500 500];
[fStack] = myCropper(I,BOX,featureFunction);
%%

%%
wholeF = '/mnt/tetra/nate/AIF/WHOLE/';
qrF = '/mnt/tetra/nate/AIF/QR/';
bioF = '/mnt/tetra/nate/AIF/BIO/';
tableF = '/mnt/tetra/nate/AIF/TAB/';

mkdir(bioF);
mkdir(qrF);
mkdir(wholeF);
mkdir(tableF);

parfor e = 1:numel(rndList)
    try
        I = imread(rndList{e});
        [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(double(I)/255,20);
        QR = imcrop(double(I)/255,qrCropBox);
        imwrite(bioStrip,[bioF num2str(e) '.tif']);
        imwrite(QR,[qrF num2str(e) '.tif']);
        imwrite(double(I)/255,[wholeF num2str(e) '.tif']);
        csvwrite([tableF num2str(e) '.csv'],qrCropBox);
    catch
        
    end
    e
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auto clicks for RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pad = 50;
QRSHEET = [];
QR_PIECES = [];

cnt = 1;
clickedFileList = {};
clickPoints = {};

for e = 882:1:numel(FileList)
    try
        I = double(imread(FileList{e}))/255;
        I(:,(end-100):end,:) = [];


        Lab = rgb2lab(I);
        RED = Lab(:,:,2) > 25.5;
        RED = imclose(RED,strel('square',51));
        RED = bwlarge(RED);
        R = regionprops(RED);
        box = R(1).BoundingBox;
        box(1:2) = box(1:2) - pad;
        box(3:4) = box(3:4) + 2*pad;
        frame = imcrop(RED,box);
        frame = imclose(frame,strel('square',11));
        holes = imfill(frame,'holes') - frame;
        smallHoles = holes - bwareaopen(holes,1000);
        frame = logical(frame + smallHoles);
        qr = imcrop(I,box);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BLUE = Lab(:,:,3) < -8;
        BLUE = imclose(BLUE,strel('square',15));
        BLUE = bwlarge(BLUE,24);
        BLUEframe = imcrop(BLUE,box);
        BLUEframe = imdilate(BLUEframe,strel('square',9));
        %BLUEframe = imdilate(BLUEframe,strel('line',13,0));
        %BLUEframe = imdilate(BLUEframe,strel('line',13,90));
        %BLUEframe = imerode(BLUEframe,strel('line',13,0));
        %BLUEframe = imerode(BLUEframe,strel('line',13,90));
        %BLUEframe = imerode(BLUEframe,strel('square',5));
        holes = imfill(BLUEframe,'holes') - BLUEframe;
        smallHoles = holes - bwareaopen(holes,1000);
        BLUEframe = logical(BLUEframe + smallHoles);
        BLUEskeleton = bwmorph(BLUEframe,'skel',inf');
        BLUEskeleton = bwmorph(BLUEskeleton,'spur',inf);
        BLUEskeleton = bwmorph(BLUEskeleton,'skel',inf');
        BLUEskeleton = imdilate(BLUEskeleton,strel('square',1));
        BLUEsquares = imfill(BLUEskeleton,'holes');
        BLUEsquares = imdilate(BLUEsquares,strel('square',4));
        dB = bwboundaries(BLUEsquares);
        boxCorners = [];
        for b = 1:numel(dB)
            K = cwtK_closed_peaks(dB{b}(1:(end-1),:),{11},50,2*10^-8);
            boxCorners(:,:,b) = flipdim(dB{b}(K.pIDX,:),2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        skeleton = bwmorph(frame,'skel',inf');
        skeleton = bwmorph(skeleton,'spur',inf);
        skeleton = bwmorph(skeleton,'skel',inf');
        branchPoints = bwmorph(skeleton,'branchpoints',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% find the point orders
        br = [];
        [br(:,2),br(:,1)] = find(branchPoints==1);
        [~,TOPPOINT_IDX] = min(br(:,2));
        [~,LEFTPOINT_IDX] = min(br(:,1));
        [~,RIGHTPOINT_IDX] = max(br(:,1));
        MIDDLEPOINT_IDX = setdiff(1:4,[TOPPOINT_IDX LEFTPOINT_IDX RIGHTPOINT_IDX]);
        redPOINTS(2,:) = br(TOPPOINT_IDX,:);
        redPOINTS(4,:) = br(LEFTPOINT_IDX,:);
        redPOINTS(5,:) = br(MIDDLEPOINT_IDX,:);
        redPOINTS(6,:) = br(RIGHTPOINT_IDX,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        skeletonFrame = skeleton;
        skeletonFrame((redPOINTS(5,2) - 50):(redPOINTS(5,2) + 50),(redPOINTS(5,1) - 50):(redPOINTS(5,1) + 50)) = 0;
        skeletonFrame_lessTOP = skeletonFrame;
        skeletonFrame_lessTOP((redPOINTS(2,2) - 50):(redPOINTS(2,2) + 50),(redPOINTS(2,1) - 50):(redPOINTS(2,1) + 50)) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trace edge from left to top
        sk = [];
        [sk(:,2),sk(:,1)] = find(skeletonFrame);
        sourcePoint = snapTo(sk,redPOINTS(4,:));
        targetPoint = snapTo(sk,redPOINTS(2,:));
        ADJ = Radjacency(sk',2);
        [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
        path = sk(pathIDX,:);
        K = cwtK_filter(path,{7});
        K.K(1:50) = 0;
        K.K(end-50:end) = 0;
        [~,ptIDX] = max(abs(K.K));
        redPOINTS(1,:) = path(ptIDX,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trace edge from top to right
        sk = [];
        [sk(:,2),sk(:,1)] = find(skeletonFrame);
        sourcePoint = snapTo(sk,redPOINTS(2,:));
        targetPoint = snapTo(sk,redPOINTS(6,:));
        ADJ = Radjacency(sk',2);
        [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
        path = sk(pathIDX,:);
        K = cwtK_filter(path,{7});
        K.K(1:50) = 0;
        K.K(end-50:end) = 0;
        [~,ptIDX] = min(K.K);
        redPOINTS(3,:) = path(ptIDX,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trace edge from left to right
        sk = [];
        [sk(:,2),sk(:,1)] = find(skeletonFrame_lessTOP);
        sourcePoint = snapTo(sk,redPOINTS(4,:));
        targetPoint = snapTo(sk,redPOINTS(6,:));
        ADJ = Radjacency(sk',2);
        [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
        path = sk(pathIDX,:);
        K = cwtK_filter(path,{7});
        K.K(1:50) = 0;
        K.K(end-50:end) = 0;
        [~,ptIDX] = min(K.K);
        f1 = imdilate(K.K,ones(100,1)) == K.K;
        f2 = K.K > .02;
        fidx = find(f1.*f2);
        redPOINTS(7:8,:) = path(fidx,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        %imshow(skeleton,[]);




        %{
        tmpQ = im2colF(qr,[21 21 3],[1 1 1]);

        if exist(qE)
            qqC = PCA_REPROJ_T(tmpQ,qE,qU);
            for p = 1:size(qqC,1)
                P(:,:,p) = col2im(qqC(p,:),[21 21],[size(qr,1) size(qr,2)]);
            end
        end

        tmpF = im2colF(double(frame),[21 21],[1 1]);
        tmpS = im2colF(double(skeleton),[21 21],[1 1]);
        kidx = find(tmpF((end-1)/2,:) == 1);
        kidxS = find(tmpS((end-1)/2,:) == 1);
        tmpQ = tmpQ(:,kidx);
        QR_PIECES = [QR_PIECES tmpQ];


        QRSHEET = [QRSHEET;imresize(qr,round([500 700]/4))];
        %}
        imshow(qr,[]);
        hold on
        plot(br(:,1),br(:,2),'bo');
        plot(redPOINTS(:,1),redPOINTS(:,2),'g*');
        CL = {'r*','m*','k*','c*'};
        for b = 1:size(boxCorners,3)
            for p = 1:size(boxCorners,1)
                plot(boxCorners(p,1,b),boxCorners(p,2,b),CL{p})
            end
        end
        blue{cnt} = boxCorners;
        clickPoints{cnt} = redPOINTS;
        clickedFileList{cnt} = imread(FileList{e});
        QR{cnt} = qr;
        cnt = cnt + 1;

        hold off;
        drawnow
    catch
        
    end
end
%% make noise images for red sparkle - ver2
close all
FilePath = qr_oPath;
QR_noiseList = {};
FileExt = {'tif'};
QR_noiseList = gdig(FilePath,QR_noiseList,FileExt,1);
SPboxSequence = [1300 1300];
SPzoomSequence = [.15];
cnt = 1;
for e = 1:numel(QR_noiseList)
    for m = 1:5
        mag = (.8*rand(1)+.6);
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = QR{e};
        
        matFile = strrep(QR_noiseList{e},'.tif','.mat');
        ob = load(matFile);
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k))/255;
        end
        %{
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        drawnow
        hold off
        %}
        NEWX(:,:,:,cnt) = nIMGx;
        %GOOD_red_corners = bsxfun(@minus,GOOD_red_corners,mean(GOOD_red_corners,1));
        NEWY(cnt,:) = GOOD_red_corners(:);%*SPzoomSequence^-1;
        cnt = cnt + 1;
    end
    e
end
%%
[qS qC qU qE qL qERR qLAM] = PCA_FIT_FULL_T(QR_PIECES,4);
%% dither QR sheets
[IND,qrMap] = rgb2ind(QRSHEET,4);
%%
qrBITS = [];
for e = 1:10
    I = double(imread(FileList{e}))/255;
    Lab = rgb2lab(I);
    RED = Lab(:,:,2) > 31;
    RED = bwlarge(RED);
    R = regionprops(RED);
    box = R(1).BoundingBox;
    box(1:2) = box(1:2) - pad;
    box(3:4) = box(3:4) + 2*pad;
    frame = imcrop(RED,box);
    qr = imcrop(I,box);
    qr = imresize(qr,[round([500 700]/4)]);
    qr = rgb2ind(qr,qrMap,'nodither');
    qr = im2col(qr,[7 7],'sliding');
    qr = squeezeQbits(qr',2,3,'uint8');
    qrBITS = [qrBITS ; qr];
end
%%
%qrBITS = qrBITS(randperm(size(qrBITS,1)),:);
%qrBITS = qrBITS(1:50000,:);
%%
close all
H1 = figure;
bitCNT = 1;
DELTA = [];
N = 10000;
STORE = {};
for y = 1:size(qrBITS,2)
    for b = 1:8
        tmp = [];
        for r = 1:20
            
            rnd = randperm(size(qrBITS,1));
            
            [X] = randBitSample(qrBITS,8,8);
            %{
            % try bllinear
            [X] = randBitSample(qrBITS,4,8);
            [X] = squeezeQbits([X X],4,3,'uint8');
            %}
            
            %{
            [X] = randBitSample(qrBITS,2,8);
            [X] = squeezeQbits([X X X X],2,3,'uint8');
            %}
            
            
            [X] = quantumEncode(8,X);
            Y = bitget(qrBITS(:,y),b);
            
            
            trainX = X(rnd(1:N),:);
            testX = X(rnd((N+1):end),:);
            

            trainY = Y(rnd(1:N));
            testY = Y(rnd((N+1):end));
            
            [STORE{y,b},DIS] = findProgram(true,2,10,trainY,trainX,10,size(X,2),3);
            
            
            
            [minDELTA,programIDX] = max(DIS(:,1));
            preY = boo(testX,STORE{y,b}(programIDX,:));
            
            testD = bitwise_hamming(preY,testY);
            
            pBits_train = minDELTA/N;
            %pBits_train = testD/(size(qrBITS,1)-N);
            pBits_test = testD/(size(qrBITS,1)-N);
            
           
            
            tmp(r,:) = [pBits_train pBits_test];
            r
        end
        
        DELTA(bitCNT,:) = mean(tmp,1);
        bitCNT = bitCNT + 1;
        figure(H1);
        plot(DELTA)
        drawnow
    end
    
end
%%
ridx = find(DELTA < .1);
qrBITS = copyPackedBits(qrBITS,ridx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clicks for RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sparkle
FilePath = '/mnt/tetra/nate/seedlingDATApile/maizeData/seedlingData/';
FileList = {};
FileExt = {'nef'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
%% cali
caliFilePath = '/mnt/tetra/nate/caliSampleRAWWHOLE/';
caliFileList = {};
FileExt = {'nef'};
caliFileList = gdig(caliFilePath,caliFileList,FileExt,1);
%% gather clicks
toc
RND = randperm(numel(FileList));
sFileList = FileList(RND(1:50));
[pointList,USE] = zoomGather_nonDB(sFileList,[100 100],true);
%% combine point list
for e = 1:numel(pointList)
    qrCenterPoint{e} = mean(pointList{e},1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clicks for RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% weave point list for zoom on corner features
featurePointLists = {};
for e = 1:numel(pointList)
    for a = 1:size(pointList{e},1)
        if e == 1
             featurePointLists{a} = [];
        end
        featurePointLists{a} = [featurePointLists{a};pointList{e}(a,:)];
    end
end
%% soom-sample features on QR sheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom for RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample for zoom sequence for whole QR
boxSequence_F = [[300 300];[75 75]];
zoomSequence_F = [.08 .5];
errorSequence_F = [50 20];
[qrX_F,qrY_F,qr_FPN] = sampleZoomSequence(sFileList,pointList,boxSequence_F,zoomSequence_F,errorSequence_F,35,[-pi/32,pi/32,100],true);
%% weave zoom samples
clear map_F beta_F WINDOW_F
for z = 1:numel(zoomSequence_F)
    for e = 1:8
        fidx = qr_FPN{z} ==e;
        featureX{z}{e} = qrX_F{z}(:,:,:,fidx);
        featureY{z}{e} = qrY_F{z}(:,fidx)';
        [map_F{e}{z},beta_F{e}{z},WINDOW_F{e}{z}] = fitbit(featureX{z}{e},featureY{z}{e},4,30);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom for RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sparkle RED CORNDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparkle whole QR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample for sparkle sequence
SPboxSequence = [1300 1300];
SPzoomSequence = [.15];
SPjiggleRad = [50];
SPjiggleNumber = 30;
SProtPara = [-pi/32,pi/32,100];
[qrSX,qrSY] = sampleSparkleJitterSequence(sFileList,pointList,qrCenterPoint,SPboxSequence,SPzoomSequence,SPjiggleRad,SPjiggleNumber,SProtPara,true);
%%
subI = [];
close all
cnt = 1;
for e = 1:numel(sFileList)
    I = double(imread(sFileList{e}))/255;
    for p = 1:8
        box = point2Box(pointList{e}(p,:),[50 50]);

        %%subI(:,:,:,e) = mimcrop(I,box,1,[]);
        subI(:,:,:,cnt) = imcrop(I,box);
        imshow(subI(:,:,:,cnt),[]);
        cnt = cnt + 1;
        drawnow
    end
end
%%
T = [];
for e = 1:size(subI,4)
    T(:,:,e) = rgb2ind(subI(:,:,:,e),nMap,'nodither');
end
sz = size(T);
T = reshape(T,[prod(sz(1:2)) sz(3)]);
nT = squeezeQbits(T,1,3,'uint8');
%% generate sparkle data for fitbit
sparkleQRY = [];
for e = 1:size(qrSY{1},3)
    tmp = qrSY{1}(:,:,e)';
    sparkleQRY(e,:) = tmp(:)';
end
%% fit bit for sparkle
%[map_S,beta_S,WINDOW_S] = fitbit(qrSX{1},sparkleQRY,2,200);
[map_S,beta_S,WINDOW_S] = fitbit(NEWX,NEWY,2,200,map_S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train sparkle - single local
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tX = [];
for e = 1:size(NEWX,4)
    %tX(:,:,1,e) = rgb2ind(imresize(qrSX{1}(:,:,:,e),.35),map_S);
    tX(:,:,1,e) = rgb2ind(imresize(NEWX(:,:,:,e),.35),map_S);
    %tX(:,:,:,e) = imresize(NEWX(:,:,:,e),.35);
end
%%
%tX = NEWX;
tY = NEWY;
%tY = sparkleQRY;
%tX = qrSX{1};
%tX = reshape(tX,[size(tX,1) size(tX,2) 1 size(tX,3)]);
layers = [imageInputLayer([size(tX,1) size(tX,2) size(tX,3)],'Normalization','None');
          convolution2dLayer([11 11],5);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(size(tY,2));
          regressionLayer();];
options = trainingOptions('sgdm',...
    'InitialLearnRate',.00001,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
%[tY,suY,ssY] = zscore(tY);
offSET = mean(tY,1);
tY = bsxfun(@minus,tY,offSET);
usparkle = trainNetwork(tX,tY,layers,options);
%% rCNN for sparkle
tX = NEWX;
tY = NEWY;
nTrain = [300 5000];
% create func
func = cFlow('gpuTrainRegression');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
% assign train data
%tY = sparkleQRY;
%tX = qrSX{1};
close all
% zscore
[tY,suY,ssY] = zscore(tY);
% number to split
N = 2100;
% split into train and test
IDX = randperm(size(tX,4));
trainX = tX(:,:,:,IDX(1:N));
trainY = tY(IDX(1:N),:);
testX = tX(:,:,:,IDX((N+1):end));
testY = tY(IDX((N+1):end),:);
% set up parameters to test with 
nHood = optimizableVariable('nHood',[7,15],'Type','integer');
nKernels = optimizableVariable('nKernels',[3,15],'Type','integer');
initLearnRate = optimizableVariable('initLearnRate',[.0001,.1]);
L2Regularization = optimizableVariable('L2Regularization',[.0001,.0008]);
Momentum = optimizableVariable('Momentum',[0,1]);
exeEnvironment = 'gpu';

para = [nHood,nKernels,initLearnRate,L2Regularization,Momentum];

[BO_sparkel,RED_corner_sparkleNET] = func(nTrain,trainX,trainY,testX,testY,exeEnvironment,para);

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
[sparkle_preY] = fitbitp(qrSX{1}(:,:,:,1),4,map_S,beta_S,WINDOW_S);
sparkle_preY = reshape(sparkle_preY,[8 2]);
imshow(qrSX{1}(:,:,:,1),[]);
hold on
plot(sparkle_preY(:,1),sparkle_preY(:,2),'g*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sparkle RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom for WHOLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample for zoom sequence for whole QR
boxSequence = [[1300 1300];[900 900]];
zoomSequence = [.08 .15];
[qrX,qrY] = sampleZoomSequence(sFileList,qrCenterPoint,boxSequence,zoomSequence,[300 100],35,[-pi/32,pi/32,100],true);
%% fit for zoom sequence for whole QR
clear map_Q beta_Q WINDOW_Q
for z = 1:numel(qrX)
    [map_Q{z},beta_Q{z},WINDOW_Q{z}] = fitbit(qrX{z},qrY{z}',2,128);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom for WHOLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull out msg
oPath = '/mnt/tetra/nate/seedlingDATApile/metaDataOutput/';
FileList = FileList(randperm(numel(FileList)));
for e = 1:numel(FileList)   
   msg =  renderMetaData_dev(FileList{e},map_Q,beta_Q,WINDOW_Q,map_S,beta_S,WINDOW_S,map_F,beta_F,WINDOW_F,oPath,false,true);
end
%%

FilePath = '/mnt/tetra/nate/seedlingDATApile/metaDataOutput/';
new_qrPileFileList = {};
FileExt = {'mat'};
new_qrPileFileList = gdig(FilePath,new_qrPileFileList,FileExt,1);
%%
close all
big = [];
h1 = figure;
h2 = figure;
dynPoints = {};
cnt = 1;
for e = 1:numel(new_qrPileFileList)
    
    q = load(new_qrPileFileList{e},'wholeSheet','sP','fidx');
    
    
    try
        
        
        basePoint = q.sP(1,:,q.fidx(end));
        dyC = mean(q.sP(:,:,q.fidx(end)),1);
        J = imread(strrep(new_qrPileFileList{e},'.mat','.tif'));
        X = rgb2ind(q.wholeSheet, map);
        mask = X == 0;
        mask = imfill(mask,'holes');
        mask = bwareaopen(mask,1000);
        mask = bwlarge(mask,24);
        mask = imclose(mask,strel('square',3));
        mask = imclose(mask,strel('line',11,0));
        mask = imclose(mask,strel('line',11,90));
        mask = imclose(mask,strel('square',7));
        mask = bwmorph(mask,'spur',inf);
       
        
        
        MK = sum(mask,2)*sum(mask,1);
        RM = regionprops(logical(MK),'PixelIdxList');
        %RM = regionprops(mask,'PixelIdxList');
        
        cPoints = detectMinEigenFeatures(mask,'FilterSize',11);
        [~,sidx] = sort(cPoints.Metric,'descend');
        cPoints = cPoints(sidx(1:(24*4)));
        
        
        
        SG = [];
        groupP = [];
        NO = false;
        for g = 1:24
            Z = zeros(size(mask));
            Z(RM(g).PixelIdxList) = 1;
            kp = [];
            for p = 1:size(cPoints.Location,1)
                kp(p) = Z(round(cPoints.Location(p,2)),round(cPoints.Location(p,1)));
            end
            SG(g) = sum(kp);
            if SG(g) == 4
                fidx = find(kp);
                subG = cPoints.Location(fidx,:);
                disG = sum(subG.*subG,2);
                % upper left corner
                [~,gidx] = min(disG);
                groupP(1,:,g) = subG(gidx,:);
                subG(gidx,:) = [];
                % upper most
                [~,gidx] = min(subG(:,2));
                groupP(2,:,g) = subG(gidx,:);
                subG(gidx,:) = [];
                % left most
                [~,gidx] = min(subG(:,1));
                groupP(3,:,g) = subG(gidx,:);
                subG(gidx,:) = [];
                groupP(4,:,g) = subG;
            else
                NO = true;
                break
            end
        end
        
        
        q.wholeSheet = q.wholeSheet / max(max(max(q.wholeSheet)));
        out = flattenMaskOverlay(double(q.wholeSheet),mask,.2,'b');
        figure(h1);
        imshow(out,[]);
        groupP = double(groupP);
        if ~NO
            hold on
            plot(cPoints.Location((1:(24*4)),1),cPoints.Location((1:(24*4)),2),'bo');
            for g = 1:size(groupP,3)
                for pt = 1:size(groupP,1)
                    plot(groupP(pt,1,g),groupP(pt,2,g),'g*')
                    text(groupP(pt,1,g)+10,groupP(pt,2,g)+10,num2str(pt))
                end
            end
            
            
            
            hold off
            drawnow
            big = [big imresize(q.wholeSheet,[400 400])];
            drawnow
        end
        
        
        
        
        if ~NO
            figure(h2)
            imshow(J,[]);
            hold on
            for g = 1:size(groupP,3)
                for pt = 1:size(groupP,1)
                    plot(groupP(pt,1,g)+basePoint(1),groupP(pt,2,g)+basePoint(2),'g*')
                    text(groupP(pt,1,g)+10+basePoint(1),groupP(pt,2,g)+10+basePoint(2),[num2str(g) '-' num2str(pt)])
                end
            end
            groupP = permute(groupP,[1 3 2]);
            dynPoints{cnt} = reshape(groupP,[4*24 2]);
            dynPoints{cnt} = bsxfun(@plus,dynPoints{cnt},basePoint);
            %dynPoints{cnt} = bsxfun(@minus,dynPoints{cnt},dyC);
            
            
            dFileList{cnt} = strrep(new_qrPileFileList{e},'.mat','.tif');
            
            dynCenterPoint{cnt} = dyC;
            cnt = cnt + 1;
        end
        
        hold off
        %break
    catch ME
        ME
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sparkle DYNAMIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample for sparkle sequence for dynamic box
DYboxSequence = [900 800];
DYzoomSequence = [.10];
DYjiggleRad = [100];
DYjiggleNumber = 30;
DYrotPara = [-pi/32,pi/32,100];
[dySX,dySY] = sampleSparkleJitterSequence(dFileList,dynPoints,dynCenterPoint,DYboxSequence,DYzoomSequence,DYjiggleRad,DYjiggleNumber,DYrotPara,true);
%% generate sparkle data for fitbit
sparkleDYY = [];
for e = 1:size(dySY{1},3)
    tmp = dySY{1}(:,:,e)';
    sparkleDYY(e,:) = tmp(:)';
end
%% fit bit for sparkle
[map_D,beta_D,WINDOW_D] = fitbit(dySX{1},sparkleDYY,2,120);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sparkle DYNAMIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom for DYNAMIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample for zoom sequence for whole QR
boxSequence_F_DY = [[90 90]];
zoomSequence_F_DY = [.5];
errorSequence_F_DY = [35];
[qrX_F_DY,qrY_F_DY,qr_FPN_DY] = sampleZoomSequence(dFileList,dynPoints,boxSequence_F_DY,zoomSequence_F_DY,errorSequence_F_DY,35,[-pi/32,pi/32,100],false);
%% weave zoom samples
clear map_F_DY beta_F)DY WINDOW_F_DY featureX_DY featureY_DY
for z = 1:numel(zoomSequence_F_DY)
    for e = 1:4
        
        %fidx = qr_FPN_DY{z} == e;
        fidx = mod(qr_FPN_DY{z},4) == (e-1);
        
        featureX_DY{z}{e} = qrX_F_DY{z}(:,:,:,fidx);
        featureY_DY{z}{e} = qrY_F_DY{z}(:,fidx)';
        %{
        for i = 1:5
            imshow(featureX_DY{z}{e}(:,:,:,i),[]);
            featureY_DY{z}{e}
        end
        %}
        
        [map_F_DY{e}{z},beta_F_DY{e}{z},WINDOW_F_DY{e}{z}] = fitbit(featureX_DY{z}{e},featureY_DY{z}{e},2,20);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom for DYNAMIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom  POTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather clicks
[pointList_POTS,USE_POTS] = zoomGather_nonDB(dFileList,[500 500],true);
%% center the pot point lists
for e = 1:numel(pointList_POTS)
    POT_CEN{e} = [];
    for p = 1:2:5
        POT_CEN{e} = [POT_CEN{e};mean(pointList_POTS{e}(p:(p+1),:))];
    end
end
%% sample for zoom sequence for POTS
boxSequencePOT = [[1300 1300];[600 600]];
zoomSequencePOT = [.08 .15];
errorPOT = [300 100];
[potX,potY] = sampleZoomSequence(dFileList,POT_CEN,boxSequencePOT,zoomSequencePOT,errorPOT,35,[-pi/32,pi/32,100],false);
%% fit for zoom sequence for POTS
clear map_POT beta_POT WINDOW_POT
for z = 1:numel(qrX)
    [map_POT{z},beta_POT{z},WINDOW_POT{z}] = fitbit(potX{z},potY{z}',2,128);
    z
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom POTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sparkle POTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stack sparkle for pots 
potFileList = {};
potCenterSample = {};
potWingSample = {};
cnt = 1;
for e = 1:numel(dFileList)
    for pot = 1:3
        idx = 2*(pot-1)+1;
        potFileList{cnt} = dFileList{e};
        potCenterSample{cnt} = POT_CEN{e}(pot,:);
        potWingSample{cnt} = pointList_POTS{e}(idx:(idx+1),:);
        cnt = cnt + 1;
    end
end
%% sample for sparkle sequence for dynamic box
POT_sparkle_boxSequence = [1500 400];
POT_sparkle_zoomSequence = [.10];
POT_sparkle_jiggleRad = [50];
POT_sparkle_jiggleNumber = 30;
POT_sparkle_rotPara = [-pi/32,pi/32,100];
[POT_S_X,POT_S_Y] = sampleSparkleJitterSequence(potFileList,potWingSample,potCenterSample,POT_sparkle_boxSequence,...
    POT_sparkle_zoomSequence,POT_sparkle_jiggleRad,POT_sparkle_jiggleNumber,POT_sparkle_rotPara,false);
%% generate sparkle data for fitbit
sparkle_POT_Y = [];
for e = 1:size(POT_S_Y{1},3)
    tmp = POT_S_Y{1}(:,:,e)';
    sparkle_POT_Y(e,:) = tmp(:)';
end
%% fit bit for sparkle
[map_SP_POT,beta_SP_POT,WINDOW_SP_POT] = fitbit(POT_S_X{1},sparkle_POT_Y,2,150);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sparkle POTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom lock feature POTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather clicks
%[pointList_CONE,USE_CONE] = zoomGather_nonDB(dFileList,[800 800],true);
%{
J = [];
for e = 1:numel(dFileList)
    tmp = imread(dFileList{e});
    J = [J tmp];
end
[X,map_CONE] = rgb2ind(J,8);
%}
%% auto click snap to cone
coneBOX = {};
for e = 1:numel(dFileList)
    J = imread(dFileList{e});
    %{
    tmp = rgb2ind(J,map_CONE);
    mask = tmp == 0;
    %}
    JH = rgb2hsv(J);
    mask = (JH(:,:,1) < .12 | JH(:,:,1) > .6) & JH(:,:,3) < .35;
    mask(1:1500,:) = 0;
    
    mask = bwlarge(mask,3);
    mask = imfill(mask,'holes');
    mask = imopen(mask,strel('square',31));
    mask = bwlarge(mask,3);
    R = regionprops(logical(mask),'PixelIdxList');
    for p = 1:3
        Z = zeros(size(mask));
        Z(R(p).PixelIdxList) = 1;
        sZ = sum(Z,2);
        gidx1 = find(sZ > 0);
        gidx1 = gidx1(1);
        gZ = gradient(sZ);
        [~,H] = max(gZ);
        
        H = H - 2;
        
        Z(1:gidx1) = 0;
        Z(H:end,:) = 0;
        
        Z = bwlarge(Z);
        Z1 = find(sum(Z,1) > 2);
        q1 = Z1(1);
        q2 = Z1(end);
        
        imshow(J,[]);
        hold on
        plot(1:size(Z,2),gidx1*ones(size(Z,2),1),'g')
        plot(1:size(Z,2),H*ones(size(Z,2),1),'g')
        plot(q1*ones(size(Z,1),1),1:size(Z,1),'g')
        plot(q2*ones(size(Z,1),1),1:size(Z,1),'g')
        
       
        
        coneBOX{e}(:,:,p) = [[q1 H];[q2 H];[q2 gidx1];[q1 gidx1]];
        pBOX = [coneBOX{e}(:,:,p);coneBOX{e}(1,:,p)];
        plot(pBOX(:,1),pBOX(:,2),'r');
        drawnow
        pause(.1)
        
        hold off
        
    end
end
%% sample for sparkle sequence for dynamic box
CONE_sparkle_boxSequence = [1500 400];
CONE_sparkle_zoomSequence = [.10];
CONE_sparkle_jiggleRad = [50];
CONE_sparkle_jiggleNumber = 30;
CONE_sparkle_rotPara = [-pi/32,pi/32,100];
[CONE_S_X,CONE_S_Y] = sampleSparkleJitterSequence(dFileList,potWingSample,potCenterSample,CONE_sparkle_boxSequence,...
    CONE_sparkle_zoomSequence,CONE_sparkle_jiggleRad,CONE_sparkle_jiggleNumber,CONE_sparkle_rotPara,true);
%% generate sparkle data for fitbit
sparkle_POT_Y = [];
for e = 1:size(POT_S_Y{1},3)
    tmp = POT_S_Y{1}(:,:,e)';
    sparkle_POT_Y(e,:) = tmp(:)';
end
%% fit bit for sparkle
[map_SP_POT,beta_SP_POT,WINDOW_SP_POT] = fitbit(POT_S_X{1},sparkle_POT_Y,2,150);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zoom lock feature POTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clip out qr code sheets
RND = 1:1:numel(FileList);
for e = 1:numel(RND)
    zoomSparkleLock_QR(e,FileList{RND(e)},boxSequence,zoomSequence,map_Q,beta_Q,WINDOW_Q,SPboxSequence,SPzoomSequence,map_S,beta_S,WINDOW_S,boxSequence_F,zoomSequence_F,map_F,beta_F,WINDOW_F)
   %zoomSparkleLock_QR();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = double(imread(FileList{3}))/255;
%I = double(imread(caliFileList{3}))/255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate grid for scatter zoom for QR
[pY,pSZ] = genIXgrid2(size(I),[800 800],[0 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first order snap to center of QR
pY = fliplr(pY);
ITER = [2 2];
[pY] = runZoomSequence(I,2,pY,boxSequence,zoomSequence,map_Q,beta_Q,WINDOW_Q,ITER);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparkle regress to red corner boxes from end points of first order snap
fP = pY(:,:,end);
[sP] = runSparkleSequence(I,2,fP,SPboxSequence,SPzoomSequence,map_S,beta_S,WINDOW_S,[8 2],usparkle,offSET);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zoom snap to red corners
ITER = [1 2];
for centerPoint = 1:size(fP,1)
    for featurePoint = 1:size(sP,1)
        [zoomTrail] = runZoomSequence(I,4,sP(featurePoint,:,centerPoint),boxSequence_F,zoomSequence_F,map_F{featurePoint},beta_F{featurePoint},WINDOW_F{featurePoint},ITER);
        sP(featurePoint,:,centerPoint) = squeeze(zoomTrail(:,:,end));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for good snaps
msg = {};
[msg] = readQRdata(I,sP,10,false);
cnt = 1;
dy_predict_points = [];
probeCP = [];
goodPoints = [];
GOOD_red_corners = [];
% find good msg
for centerPoint = 1:numel(msg)
    if ~isempty(msg{centerPoint})
        goodPoints = [goodPoints;centerPoint];
        probeCP(cnt,:) = fP(centerPoint,:);
        GOOD_red_corners = cat(3,GOOD_red_corners,sP(:,:,centerPoint));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sparkle predict the blue dynamic boxes
        %[dy_predict_points(:,:,cnt)] = runSparkleSequence(I,2,probeCP(cnt,:),DYboxSequence,DYzoomSequence,map_D,beta_D,WINDOW_D,[24*4 2]);
        cnt = cnt + 1;
    end
end
%% make noise images for red sparkle
close all
FilePath = qr_oPath;
QR_noiseList = {};
FileExt = {'tif'};
QR_noiseList = gdig(FilePath,QR_noiseList,FileExt,1);
SPboxSequence = [1300 1300];
SPzoomSequence = [.15];
cnt = 1;
for e = 1:numel(QR_noiseList)
    for m = 1:5
        mag = (.8*rand(1)+.6);
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = imread(QR_noiseList{e});
        matFile = strrep(QR_noiseList{e},'.tif','.mat');
        ob = load(matFile);
        
        
        nIMGx = rand([fliplr(SPboxSequence) 3]);
        nIMGx = imresize(nIMGx,SPzoomSequence);
        GOOD_red_corners = bsxfun(@minus,ob.GOOD_red_corners(:,:,end),ob.GOOD_red_corners(1,:,end));
        GOOD_red_corners = SPzoomSequence_tmp*bsxfun(@plus,GOOD_red_corners,[26 26]);
        
        
        
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k))/255;
        end
        %{
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        drawnow
        hold off
        %}
        NEWX(:,:,:,cnt) = nIMGx;
        %GOOD_red_corners = bsxfun(@minus,GOOD_red_corners,mean(GOOD_red_corners,1));
        NEWY(cnt,:) = GOOD_red_corners(:);%*SPzoomSequence^-1;
        cnt = cnt + 1;
    end
    e
end
%%
POTP = dy_predict_points;
dy_predict_points = mean(dy_predict_points,3);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zoom lock on blue boxes
ITER = [4];
BLUEzoomTrail = [];
for b = 1:size(dy_predict_points,1)
    level = mod(b,4)+1;
    [BLUEzoomTrail(b,:,:)] = runZoomSequence(I,2,dy_predict_points(b,:),boxSequence_F_DY,zoomSequence_F_DY,map_F_DY{level},beta_F_DY{level},WINDOW_F_DY{level},ITER);
    dy_predict_points(b,:) = squeeze(BLUEzoomTrail(b,:,end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate grid for scatter zoom for pots
[potY,potSZ] = genIXgrid2(size(I),[400 400],[0 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first order snap to POTS
potY = fliplr(potY);
ITER = [2 2];
[potY] = runZoomSequence(I,2,potY,boxSequencePOT,zoomSequencePOT,map_POT,beta_POT,WINDOW_POT,ITER);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparkle regress to red corner boxes from end points of first order snap
[sparklePOT_YLINE] = runSparkleSequence(I,2,potY(:,:,end),POT_sparkle_boxSequence,POT_sparkle_zoomSequence,...
    map_SP_POT,beta_SP_POT,WINDOW_SP_POT,[2 2]);
%% view testers
close all
imshow(I,[]);
hold on
% draw paths for zoom sequence for QR sheet
plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'g*')
plot(squeeze(pY(:,1,end)),squeeze(pY(:,2,end)),'r*')
for p = 1:size(pY)
    plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k')
    plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k.')
    hold on
end
% draw web for red corners
for p = 1:size(pY)
    plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k')
    plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k.')
    hold on
    plot(squeeze(pY(p,1,end)),squeeze(pY(p,2,end)),'r*')
    plot(squeeze(pY(p,1,end)),squeeze(pY(p,2,end)),'ko')

    plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'g*')
    plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'ko')
end
for p = 1:size(sP,3)
    plot(sP(:,1,p),sP(:,2,p),'k.')
    plot(sP(:,1,p),sP(:,2,p),'ko')
    VV = repmat(squeeze(pY(p,:,end)),[size(sP) 1]);
    for q = 1:size(sP,1)
        plot([VV(q,1);sP(q,1,p)],[VV(q,2);sP(q,2,p)],'r')
    end
end
%%
% draw paths for zoom sequence for QR sheet
plot(squeeze(potY(:,1,1)),squeeze(potY(:,2,1)),'g*')
plot(squeeze(potY(:,1,end)),squeeze(potY(:,2,end)),'r*')
%%
for p = 1:size(potY)
    plot(squeeze(potY(p,1,:)),squeeze(potY(p,2,:)),'m')
    plot(squeeze(potY(p,1,:)),squeeze(potY(p,2,:)),'m.')
    hold on
end
for p = 1:size(sparklePOT_YLINE,3)
    plot(sparklePOT_YLINE(:,1,p),sparklePOT_YLINE(:,2,p),'m')
    plot(sparklePOT_YLINE(1,1,p),sparklePOT_YLINE(1,2,p),'m*')
    plot(sparklePOT_YLINE(1,1,p),sparklePOT_YLINE(1,2,p),'ko')
    plot(sparklePOT_YLINE(2,1,p),sparklePOT_YLINE(2,2,p),'m*')
    plot(sparklePOT_YLINE(2,1,p),sparklePOT_YLINE(2,2,p),'ko')
end


%%
for p = 1:size(sP,3)
    
    plot(sP(:,1,p),sP(:,2,p),'k.')
    
    plot(sP(:,1,p),sP(:,2,p),'k.')
    
end
plot(probeCP(:,1),probeCP(:,2),'b*')
for e = 1:size(dy_predict_points,3)
    plot(dy_predict_points(:,1),dy_predict_points(:,2),'b.')
    plot(dy_predict_points(:,1),dy_predict_points(:,2),'ro')
end
for e = 1:size(BLUEzoomTrail,1)
    plot(squeeze(BLUEzoomTrail(e,1,:)),squeeze(BLUEzoomTrail(e,2,:)),'b')
end
plot(squeeze(BLUEzoomTrail(:,1,1)),squeeze(BLUEzoomTrail(:,2,1)),'b.')
plot(squeeze(BLUEzoomTrail(:,1,1)),squeeze(BLUEzoomTrail(:,2,1)),'go')
%% for pretty
close all
imshow(I,[]);
hold on
squareQR = sP([1 2 5 4 1],:,goodPoints(1));
squareHU = sP([2 3 6 5 2],:,goodPoints(1));
squareDY = sP([4 6 8 7 4],:,goodPoints(1));
%%
plot(squareQR(:,1),squareQR(:,2),'k','LineWidth',1)
%%
plot(squareHU(:,1),squareHU(:,2),'g','LineWidth',1)
plot(squareDY(:,1),squareDY(:,2),'b','LineWidth',1)
%%
hold on
plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'g*')
plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'ko')
%%
p = 15
plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k')
plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k.')
hold on
plot(squeeze(pY(p,1,end)),squeeze(pY(p,2,end)),'r*')
plot(squeeze(pY(p,1,end)),squeeze(pY(p,2,end)),'ko')
%%
for p = 1:24
plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k')
plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k.')
hold on
plot(squeeze(pY(p,1,end)),squeeze(pY(p,2,end)),'r*')
plot(squeeze(pY(p,1,end)),squeeze(pY(p,2,end)),'ko')

plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'g*')
plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'ko')
end
%%

%for p = 1:size(sP,3)
p = 6;
plot(sP(:,1,p),sP(:,2,p),'k.')
plot(sP(:,1,p),sP(:,2,p),'ko')
VV = repmat(squeeze(pY(p,:,end)),[size(sP) 1]);
for q = 1:size(sP,1)
    plot([VV(q,1);sP(q,1,p)],[VV(q,2);sP(q,2,p)],'k')
end
    %%
    
plot(probeCP(:,1),probeCP(:,2),'b*')
for e = 1:size(dy_predict_points,3)
    plot(dy_predict_points(:,1),dy_predict_points(:,2),'b.')
    plot(dy_predict_points(:,1),dy_predict_points(:,2),'ko')
end
%%
%for p = 1:size(sP,3)

    
%end
%%
plot(squeeze(pY(:,1,end)),squeeze(pY(:,2,end)),'r*')
for p = 1:size(pY)
    plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k')
    plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'b.')
    hold on
end
for p = 1:size(sP,3)
    
    plot(sP(:,1,p),sP(:,2,p),'k.')
    
    plot(sP(:,1,p),sP(:,2,p),'k.')
    
end







plot(probeCP(:,1),probeCP(:,2),'b*')
for e = 1:size(dy_predict_points,3)
    plot(dy_predict_points(:,1),dy_predict_points(:,2),'b.')
end
for e = 1:size(BLUEzoomTrail,1)
    plot(squeeze(BLUEzoomTrail(e,1,:)),squeeze(BLUEzoomTrail(e,2,:)),'b')
end
