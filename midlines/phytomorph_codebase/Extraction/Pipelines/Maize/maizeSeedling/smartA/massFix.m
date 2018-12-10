clear all
%%
FilePath{1} = '/mnt/tetra/nate/MN_MAT_LIST/All_030218-2018-03-02-20-00-22.9/';
FilePath{2} = '/mnt/tetra/nate/caliSample/';
for fl = 1:numel(FilePath)
    FileList{fl} = {};
    FileExt = {'mat'};
    FileList{fl} = gdig(FilePath{fl},FileList{fl},FileExt,1);
end
%% gather NEF files
imageFilePath = {};
imageFilePath{1} = '/mnt/tetra/nate/hirschSampleRAWWHOLE/';
imageFilePath{2} = '/mnt/tetra/nate/caliSampleRAWWHOLE/';
for fl = 1:numel(FilePath)
    imageFileList{fl} = {};
    imageFileExt = {'nef'};
    imageFileList{fl} = gdig(imageFilePath{fl},imageFileList{fl},imageFileExt,1);
    imageFileList{fl} = imageFileList{fl}(randperm(numel(imageFileList{fl})));
end
%% make whole list
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
N = 500;
%totD = 100;
n = 1;
parfor fl = 1:numel(FileList)
    totD = numel(FileList{fl});
    totD = 500;
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
%%
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
%%
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
%% try it on crop boxes
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
    























