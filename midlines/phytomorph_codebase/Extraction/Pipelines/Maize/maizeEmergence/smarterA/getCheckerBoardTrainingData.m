%% try guosheng
FilePath = '/mnt/spaldingdata/Red_cap_popper/Images_Testing/IMG_0273/'; % work
%FilePath = '/mnt/spaldingdata/Red_cap_popper/Images_Testing/IMG_0282/'; % fail
%FilePath = '/mnt/spaldingdata/Red_cap_popper/Images_Testing/IMG_0276/'; % fail
%FilePath = '/mnt/spaldingdata/Red_cap_popper/Images_Testing/IMG_0280/'; % fail

FileList = {};
FileExt = {'CR2'};
FileList = gdig(FilePath,FileList,FileExt,1);
[tform] = getCameraPara(FileList{1});
CB = getRectifiedImage(FileList{2},tform);
[boundingBox,centerPoints,MASK,I] = getCropBoxes(FileList,tform,168,10,150,1);
[boundingBox,centerPoints,LABELS] = orderCropBoxes(MASK,boundingBox,centerPoints,I,true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct file list for training
% 1: get list of all files from user = jgustin
%       - order via sdig style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: gather data from CyVerse
% this is the RILs only path
dataPath = '/iplant/home/jgustin%'; 
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
%% step 2: parse data
FileList = {};
FileExt = {'tiff'};
fprintf(['Found:' num2str(numel(r)) ' records \n']);
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%% step 3: make stacks 
[FileList] = orderFrom_gdig(FileList,{});
%% step 4: sort stacks
rm = [];
for e = 1:numel(FileList)
    nm = [];
    n = [];
    try
        for s = 1:numel(FileList{e})
            [~,nm] = fileparts(FileList{e}{s});
            n(s) = str2num(nm);
        end

        [~,sidx]= sort(n);
        FileList{e} = FileList{e}(sidx);
    catch
        rm = [rm e];
    end
end
FileList(rm) = [];
%%
BKFileList = FileList;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read master list and intersect with file list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read master list measurements for emergence
ml = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/EmergenceAssay_Master_list.csv');
toScan = {};
iPlantDirectory = 29;
for e = 2:size(ml,1)
    if ~isempty(ml{e,iPlantDirectory})
        if ischar(ml{e,iPlantDirectory})
            toScan{end+1} = ml{e,iPlantDirectory};
        else
            ml{e,iPlantDirectory}
            e
        end
    end
end
toScan = unique(toScan);
%% get intersection with masterlist data
CMD = 'iget -f /iplant/home/jgustin/maizeData/coleoptileEmergence/handScores/handscoreTrainSet.csv /home/nate/Downloads';
system(CMD)
D = readtext('/home/nate/Downloads/handscoreTrainSet.csv');
scoreNUMF = 250;
import java.util.Map.*;
import java.util.HashMap;
handData_t0 = HashMap();
handData_t1 = HashMap();
handData_t2 = HashMap();
SCORED = {};
for e = 2:size(D,1)
    emerFr = D{e,5};
    gidx = strfind(D{e,1},filesep);
    key = [lower(D{e,1}((gidx(end)+1):(end))) '-' lower(D{e,2}(1)) '-' lower(D{e,3}(1)) lower(num2str(D{e,4}))];  
    SCORED{end+1} = D{e,1};

    handData_t0.put(key,emerFr);
    eF = zeros(scoreNUMF,1);
    if ~isnan(emerFr) & ~isinf(emerFr) & emerFr <= 250
        eF(emerFr) = 1;
    end
    handData_t1.put(key,eF);
    handData_t2.put(key,cumsum(eF));
end
%% get intersection with masterlist data - OLD WAY
rootP = {};
for e = 1:numel(FileList)
    [rootP{e}] = fileparts(FileList{e}{1});
end
[U,sidx,~] = intersect(rootP,toScan);
FileList = FileList(sidx);
%{
%% load hand score
scoreNUMF = 250;
em = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/Emergence_hand_score_merged.csv');
iPlantDirectory = 29;
iplantName = 27+1;
iplantName = 29;
hcName = 20+1+1;
posName = 10+1+1;
genoName = 7+1+1;
FRAME_SCORE = 34;
SCORED = {};

import java.util.Map.*;
import java.util.HashMap;
handData_t0 = HashMap();
handData_t1 = HashMap();
handData_t2 = HashMap();
for e = 2:size(em,1)
    e
    gidx = strfind(em{e,iplantName},filesep);
    key = [lower(em{e,iplantName}((gidx(end)+1):(end))) '-' lower(em{e,hcName}(1)) '-' lower(em{e,29+1}(1)) lower(num2str(em{e,30+1}))];  
    SCORED{end+1} = em{e,iplantName};
    fidx = strfind(SCORED{end},filesep);
    key = lower([SCORED{end}((fidx(end)+1):end) '-' em{e,21+1}(1) '-' em{e,29+1} num2str(em{e,30+1})]);
    
    emerFr = em{e,end};
    %{
    if strcmp(lower(['20170220_Camera1-' LABELS{22}]),key)
        emerFr = 71;
    end
    if strcmp(lower(['20170220_Camera1-' LABELS{50}]),key)
        emerFr = 69;
    end
    if strcmp(lower(['20170220_Camera1-' LABELS{12}]),key)
        emerFr = 85;
    end
    if strcmp(lower(['20170220_Camera1-' LABELS{131}]),key)
        emerFr = 201;
    end
    
    
    if strcmp(lower(['20170220_Camera2-' LABELS{111}]),key)
        emerFr = 236;
    end
    if strcmp(lower(['20170220_Camera2-' LABELS{15}]),key)
        emerFr = 117;
    end
    if strcmp(lower(['20170220_Camera2-' LABELS{110}]),key)
        emerFr = 76;
    end
    if strcmp(lower(['20170220_Camera2-' LABELS{4}]),key)
        emerFr = 98;
    end
   %}
    
    
    if ischar(emerFr)
        emerFr = str2num(emerFr);
    end
    handData_t0.put(key,emerFr);
    eF = zeros(scoreNUMF,1);
    if ~isnan(emerFr) & ~isinf(emerFr) & emerFr <= 250
        eF(emerFr) = 1;
    end
    handData_t1.put(key,eF);
    handData_t2.put(key,cumsum(eF));
    
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intersect the file list with  the scored list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UQ = unique(SCORED);
rootP = {};
for e = 1:numel(FileList)
    [rootP{e}] = fileparts(FileList{e}{1});
end
[U,sidx,~] = intersect(rootP,UQ);
FileList = FileList(sidx);
%% look for matches in file list- SEARCH
%toMatch= {'20170516_Camera3','20170608_Camera2','20170613_Camera4'};
toMatch1 = {'20170613_Camera2','20170613_Camera3'};
toMatch1 = {'20170131_Camera3'};
toMatch1 = {'20170613_Camera2'};
kidx = zeros(size(FileList,1),1);
for e = 1:numel(FileList)
    for m = 1:numel(toMatch1)
        if ~isempty(strfind(FileList{e}{1},toMatch1{m}))
            kidx(e) = 1;
        end
    end
end
%FileList = FileList(find(kidx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract to disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesToMeasure = 250;
oPath = '/mnt/tetra/nate/forEM2';
for e = 1:2
       
        [tform] = getCameraPara(FileList{e}{1});
        
        
        CB = getRectifiedImage(FileList{e}{2},tform);
        
        %{
        I = imread(FileList{e}{1});
        imshow(I,[]);
        hold on
        for p = 1:size(imagePoints,1)
            text(imagePoints(p,1),imagePoints(p,2),num2str(p),'BackGround','r')
        end
        drawnow
        %}
        
        [boundingBox,centerPoints,MASK,I] = getCropBoxes(FileList{e},tform,168,10,150,1);
        [boundingBox,centerPoints,LABELS] = orderCropBoxes(MASK,boundingBox,centerPoints,I,true);
        cropTimeSeriesFromRectified(FileList{e},boundingBox,[1:1:framesToMeasure],LABELS,tform,oPath,[],[],101,[],[],[],[],[],[],[],[],[],[],[],[]);
%cropTimeSeriesFromRectified(FileList,boundingBox,timeSeries,LABELS,tform,storePath,emergenceNet,BESTemergence,nsz,zU,zE,zU2,zE2,WINDOW_SZ,Zmean,Zstd,CELLMASK,NUMPC,nI,subBOX,Tscore)
end
%% 
%% BIG NEXT
FilePath = '/mnt/tetra/nate/forEM2/';
tifFileList = {};
FileExt = {'tif','TIF'};
tifFileList = gdig(FilePath,tifFileList,FileExt,1);
%% split into ordered sets by key
imageKey = {};
for e = 1:numel(tifFileList)
    [p,n,ext] = fileparts(tifFileList{e});
    fidx = strfind(n,'--');
    imageKey{e} = n(1:(fidx(1)-1));
end
uniqueImageKey = unique(imageKey);
%{
%% check if any do not have 250 frames
for e = 1:numel(uniqueImageKey)
    parfor f = 1:numel(tifFileList)
        if ~isempty(strfind(tifFileList{f},uniqueImageKey{e}))
            grp(f) = e;
        end
    end
    e
end
nUQ = unique(grp);
rm = [];
for e = 1:numel(nUQ)
    if sum(grp==nUQ(e)) ~= 249
        rm = [rm e];
    end
end
%%
for e = 1:numel(rm)
    fidx = find(grp==rm(e));
    for f = 1:numel(fidx)
        delete(tifFileList{fidx(f)})
    end
end
%}
%% loader from disk
nsz = 151;



TYPE = zeros(numel(tifFileList),1);
TYPEBEST = TYPE;
VALUE = zeros(numel(tifFileList)/249,1);
cnt = 1;
cnt2 = 1;
cnt3 = 1;
CELLMASK = zeros(nsz,nsz);
CELLMASK((nsz-1)/2,(nsz-1)/2) = 1;
CELLMASK = bwdist(CELLMASK);
CELLMASK = double(CELLMASK < .9*(nsz-1)/2);
RC = regionprops(logical(CELLMASK),'BoundingBox');
BOX = RC(1).BoundingBox;
junk = imcrop(CELLMASK,BOX);
tiffLOAD = 1;
if tiffLOAD
    rZ = zeros([size(junk,1) size(junk,2) 3 numel(tifFileList)]);
end
TYPEF = [];
for e = 1:numel(uniqueImageKey)
    
    baseName = uniqueImageKey{e};
    if tiffLOAD
        parfor fr = 1:249
            fileName = [FilePath baseName '--f' num2str(fr) '.tif'];
            Z{fr} = imresize(double(imread(fileName)),[nsz nsz]);
            Z{fr} = bsxfun(@times,Z{fr},CELLMASK);
            Z{fr} = imcrop(Z{fr},BOX);
            %nm
        end
    end
    
    
    YV = handData_t2.get(uniqueImageKey{e});
    BB = handData_t1.get(uniqueImageKey{e});
    
    BBC = imdilate(BB,strel('disk',5));
    BBC2 = imdilate(BB,strel('disk',9));
    %{
    R = regionprops(logical(BBC2-BBC),'PixelIdxList');
    BB = zeros(size(BBC));
    if ~isempty(R)

        BB1 = zeros(size(BBC));
        BB2 = zeros(size(BBC));
        BB1(R(1).PixelIdxList) = 1;
        if numel(R) == 2
            BB2(R(2).PixelIdxList) = 1;
        end
        BB = BB1 + BBC*2 + 3*BB2;
    end
    %}
    BB = BBC2;
    Yv = handData_t0.get(uniqueImageKey{e});
    for fr = 1:249
        TYPE(cnt2) = YV(fr);
        cnt2 = cnt2 + 1;
    end
    
    TYPEF = [TYPEF;YV'];
    
    for fr = 1:249
        TYPEBEST(cnt3) = BB(fr);
        cnt3 = cnt3 + 1;
    end
    
    
    
    %TYPE = [TYPE ; YV];
    VALUE(e) = Yv;
    
    if tiffLOAD
        for fr = 1:249
            rZ(:,:,:,cnt) = Z{fr};
            cnt = cnt + 1;
        end
    end
    
    
    e
end
%% make histograms for Value
parfor e = 1:size(rZ,4)
    tmp = rgb2hsv_fast(rZ(:,:,:,e)/255,'','V');
    [H(:,e)] = hist(tmp(:),linspace(0,1,256));
    e
end
%% select histogram from first stack first cell
uH = mean(H(:,1:249),2);
%% correct images
%crZ = rZ;
CELLMASK_sub = imcrop(CELLMASK,BOX);
cidx = find(CELLMASK_sub);


toMatch = rgb2hsv_fast(mean(rZ(:,:,:,1:10:end)/255,4),'','V');
  
  
parfor e = 1:size(rZ,4)
    tmp = rgb2hsv_fast(rZ(:,:,:,e)/255);
    
    
    tmp(:,:,3) = imhistmatch(tmp(:,:,3),toMatch);
    
    
    rZ(:,:,:,e) = 255*hsv2rgb(tmp);
    e
    % imshow(cat(2,rZ(:,:,:,e)/255,tmp),[]);
end
%% subtract first 20 frames
rZtmp = reshape(rZ,[size(rZ,1) size(rZ,2) 3 249 size(rZ,4)/249]);
for e = 1:size(rZtmp,5)
    tmpU = mean(rZtmp(:,:,:,1:20,e),4);
    rZtmp(:,:,:,:,e) = bsxfun(@minus,rZtmp(:,:,:,:,e),tmpU);
    e
end
rZtmp = reshape(rZtmp,size(rZ));
%%
close all
SET = 50;
for e = 1:249
    img = rZtmp(:,:,:,e,SET)/255;
    imshow(img,[0 1])
    drawnow
end
%% view one of the corrected images
close all
imshow(rZ(:,:,:,15000)/255,[])
%% decompse loaded data
sz = size(rZ);
nZ = reshape(rZ,[prod(sz(1:3)) prod(sz(4))]);
[imgSIM,zC,zU,zE] = PCA_FIT_FULL_T(nZ,25);
%% reshape loaded data
%zCfinal = zC;
SELDIMS = 1:5;
[zCfinal Zmean Zstd] = zscore(zC(SELDIMS,:),1,2);
%zCfinal = [zCfinal;zC];
Hsz = size(zCfinal);
zCfinal = reshape(zCfinal,[Hsz(1) 249 size(zCfinal,2)/249]);
zCfinal = reshape(zCfinal,[size(zCfinal,1) size(zCfinal,2) 1 size(zCfinal,3)]);


rCfinal = reshape(zC(SELDIMS,:),[Hsz(1) 249 size(zC,2)/249]);
rCfinal = reshape(rCfinal,[size(rCfinal,1) size(rCfinal,2) 1 size(rCfinal,3)]);
%% setup rClass
BEST = [];
WINDOW_SZ = 11;
str = 1;
HUM = [];
for tr = 1:size(zCfinal,4)
    stp = str + 249-1;
    tmpS = im2col(zCfinal(:,:,1,tr),[size(zCfinal,1) WINDOW_SZ]);
    tmpS = reshape(tmpS,[size(zCfinal,1) WINDOW_SZ size(tmpS,2)]);
    
    tmpSz = im2col(rCfinal(:,:,1,tr),[size(rCfinal,1) WINDOW_SZ]);
    tmpSz = reshape(tmpSz,[size(rCfinal,1) WINDOW_SZ size(tmpSz,2)]);
    
    Ze = zeros(size(tmpS));
    
    tmpS = reshape(tmpS,[size(tmpS,1) size(tmpS,2) 1 size(tmpS,3)]);
    tmpSz = reshape(tmpSz,[size(tmpSz,1) size(tmpSz,2) 1 size(tmpSz,3)]);
    Ze = reshape(Ze,[size(Ze,1) size(Ze,2) 1 size(Ze,3)]);
    GG = cat(3,tmpS,tmpSz,Ze);
    
    BEST = cat(4,BEST,GG);
    
    tmpBB = TYPEBEST(str:stp);
    tmpBB = tmpBB((WINDOW_SZ-1)/2:(end-(WINDOW_SZ-1)/2-1));
    
    HUM = [HUM;tmpBB];
    tr
    str = stp + 1;
end
%BEST = reshape(BEST,[size(BEST,1) size(BEST,2) 1 size(BEST,3)]);
%% train dynamic
layers = [ ...
    imageInputLayer([size(BEST,1) size(BEST,2) 3])
    convolution2dLayer([1 11],15)
    reluLayer()
    %maxPooling2dLayer([1 3],'Stride',1);
    %convolution2dLayer([1 3],15)
    %reluLayer()
    convolution2dLayer([size(BEST,1) 1],5)
    reluLayer()
    %maxPooling2dLayer([3 1])
    fullyConnectedLayer(2)
    softmaxLayer()
    classificationLayer()];
options = trainingOptions('sgdm','MaxEpochs',8,'InitialLearnRate',0.00005,'ExecutionEnvironment','parallel','Plots','training-progress');
emergenceDynamic = trainNetwork(single(BEST),categorical(HUM),layers,options);
%% my metrics implement - frenet
K = [];
for e = 1:size(rCfinal,4)
    K(:,:,:,e) = myFrenet(rCfinal(:,:,1,e),7,@(X)funcG(X));
    e
end
K = abs(K);
ksz = size(K);
K = reshape(K,[ksz(1) prod(ksz(2:4))]);
[K,Kmu,Ksigma] = zscore(K,1,2);
K = reshape(K,ksz);
%% LSTM
layers = [ ...
    sequenceInputLayer(5)
    lstmLayer(5)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
options = trainingOptions('sgdm','MaxEpochs',600,'InitialLearnRate',0.0005,'Plots','training-progress');
toA = 20;
for e = 1:size(zCfinal,4)
%{
    nX{e} = [single(zCfinal(:,:,:,e));single(rCfinal(:,:,:,e))];
    init = mean(nX{e}(:,toA),2);
    nX{e} = bsxfun(@minus,nX{e},init);
    imfilter(nX{e},fspecial('average',[1 5]),'replicate');
%}
    nX{e} = K(:,:,:,e);
end
for e = 1:(size(TYPEF,1))
    YY{e} = categorical(TYPEF(e,2:end));
end
%%
mag = max(K(:));
for e = 1:numel(nX)
    plot(nX{e}');
    hold on
    axis([0 250 -mag mag])
    plot(mag*(double(YY{e})-1),'r')
    hold off
    drawnow
    pause(.2)
    
end
%%
emergenceDynamic = trainNetwork(nX,YY,layers,options);

%%
optimVars = [
    optimizableVariable('lstmStates',[2 20],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
maxE = 10;
toSlow = 1;
exeEnvironment = {'cpu',true};


optVarsTest.lstmStates = 2;
optVarsTest.InitialLearnRate = 1e-3;
optVarsTest.Momentum = .8;
optVarsTest.L2Regularization = 1e-10;



SLAP = 750;
[valError,cons] = makeObjFcn_emerge(...
                                    optVarsTest,...
                                    nX(1:SLAP),...
                                    YY(1:SLAP),...
                                    nX((SLAP+1):end),...
                                    YY((SLAP+1):end),...
                                    toSlow,...
                                    'none',...
                                    maxE,...
                                    'cpu');

maxE = 400;
toSlow = 1;
exeEnvironment = {'gpu',false};
[BayesObject] = hyperPdeploy_emerge(optimVars,nX(1:SLAP),...
                           YY(1:SLAP),...
                            nX((SLAP+1):end),...
                            YY((SLAP+1):end),...
                            exeEnvironment,maxE,toSlow);
%% remote GPU - train on GPU via condor and hyper parameters
func = cFlow('hyperPdeploy_emerge');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
% max function evaluations
maxE = 200;
% slow down training to get stable error
toSlow = 1;
% 
maxEval = 10000;
% max time for training hyper parameters
maxTime = 2*60*60;
% execution environment gpu
exeEnvironment = {'gpu',false};
% setup the variables to optimize
optimVars = [
    optimizableVariable('lstmStates',[2 20],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
% call the func and store the result in beta0 = b0
bO = func(optimVars,nX(1:SLAP),...
                           YY(1:SLAP),...
                            nX((SLAP+1):end),...
                            YY((SLAP+1):end),...
                            exeEnvironment,maxE,toSlow,maxE,maxTime);

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
hPara = cFlowLoader(bO);
%%
exeEnvironment = {'cpu',true};


optimVars = [
    optimizableVariable('lstmStates',[2 20],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];

bO_local = hyperPdeploy_emerge(optimVars,nX(1:SLAP),...
                           YY(1:SLAP),...
                            nX((SLAP+1):end),...
                            YY((SLAP+1):end),...
                            exeEnvironment,maxE,toSlow,maxE,maxTime);


%%
[valError,cons,trainedNet] = makeObjFcn_emerge(...
                                    hPara.XAtMinEstimatedObjective,...
                                    nX,...
                                    YY,...
                                    '',...
                                    '',...
                                    0,...
                                    'training-progress',...
                                    200,...
                                    'cpu');
%%
func = cFlow('myGPUtest');
func.setMCRversion('v930');
func.setGPU(1);

layers = [ ...
    sequenceInputLayer(5)
    lstmLayer(5)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

net = func(nX,YY,layers,'gpu',500);


auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%%

%% stack for hmm

EXTRA = predict(emergenceDynamic,single(BEST));
hmmBEST = BEST;
hmmBEST(:,:,3,:) = [];
hmmBEST = permute(hmmBEST,[1 3 2 4]);
Bsz = size(hmmBEST);

hmmBEST = reshape(hmmBEST,[prod(Bsz(1:2))  Bsz(3) Bsz(4)]);
fidx0 = find(HUM==0);
fidx1 = find(HUM==1);
[hmm] = makeChainRepeat_foEmergence(U0,C0,U1,C1,D,HARDLINE_HOLD,n,nCOMP,gmmNUM);

%% condor apply - local
for e = [13:numel(FileList)]
    try
       [finalScore2{e}] = detectEmergence(FileList{e},'',trainedNet,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Kmu,Ksigma,Tscore); 
        
    catch ME
        
    end
end
%%
%{
%% backgroud adjustment
parfor e = 1:numel(FileList{8})
    I{e} = imread(FileList{8}{e});
    e
end
%% crop out small wood background
[sub BOX] = imcrop(I{1});
%%
close all
subI = [];
for e = 1:numel(I)
    subI(:,:,:,e) = imcrop(I{e},BOX);
    imshow(subI(:,:,:,e)/255,[]);
    title(num2str(e))
    drawnow
end
%%
close all
imshow(I{2},[])
%%
[optimizer, metric] = imregconfig('monomodal');
parfor e = 3:size(subI,4)
    tform{e} = imregtform(rgb2gray(subI(:,:,:,e)/255), rgb2gray(subI(:,:,:,2)/255), 'translation', optimizer, metric);
    e
end
%%
close all
N = 200;
outI = [];
sz = size(subI);

for k = 1:3
    outI(:,:,k) = imwarp(subI(:,:,k,N),tform{N},'OutputView',imref2d(sz(1:2)));
end
imshow(cat(3,subI(:,:,1:2,2),outI(:,:,3))/255)
%%
close all
sz = size(I{2});
for N = 3:numel(I)
    for k = 1:3
        nI{N}(:,:,k) = imwarp(I{N}(:,:,k),tform{N},'OutputView',imref2d(sz(1:2)));
    end
    N
end
imshow(cat(3,I{2}(:,:,1:2),nI{200}(:,:,3)))
%%
close all
imshow(cat(3,subI(:,:,1:2,2),subI(:,:,3,end)))
%}
%% deploy

cornPopper = @(X)detectEmergence(X,'',trainedNet,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Kmu,Ksigma,Tscore);
pF = partialFunction(cornPopper,'cornPopperNNapp');
pF.publish();
%%

%% FINAL GRADE
DELTA = {};
for e = 1:13%numel(finalScore2)
    if ~isempty(finalScore2{e})
        Tnum = e;
        [Tpth] = fileparts(FileList{Tnum}{1});
        fidx = strfind(Tpth,filesep);
        % make labels
        LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
        NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
        HL = {'d' 'p'};
        LABELS = {};
        for e1 = 1:numel(HL)
            for e2 = 1:numel(LL)
                for e3 = 1:numel(NL)
                    LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
                end
            end
        end
        
        
         TkeyT = lower([Tpth((fidx(end)+1):end) '-' LABELS{1}]);
         sc = handData_t2.get(TkeyT);
        if ~isempty(sc)
            LABELS = lower(LABELS);
            TESTValues = [];
            for l = 1:numel(LABELS)
                Tkey{l} = lower([Tpth((fidx(end)+1):end) '-' LABELS{l}]);
                sc = handData_t2.get(Tkey{l});

                sc = sc(1:min(250,numel(FileList{e})));
                TESTValues = [TESTValues ; sc];
            end

            TESTValues = reshape(TESTValues,[min(250,numel(FileList{e})) 168]);
            Tscore = [];
            for t = 1:size(TESTValues,2)
                fidx = find(TESTValues(:,t));
                if ~isempty(fidx)
                    Tscore(t) = fidx(1);
                else
                    Tscore(t) = 0;
                end
            end
            GGG{e} = Tscore;
            DELTA{e} = Tscore - finalScore2{e};
        end
        
    end
end
%% stats on DELTA
MT = [];
P = [];
close all

for e = 1:numel(DELTA)
    if ~isempty(DELTA{e})
        tmp = DELTA{e};
        tmp(finalScore2{e}==0) = 0;
        tmp2 = GGG{e};
        tmp2(finalScore2{e}==0) = 0;
        tmp3 = finalScore2{e};
        tmp3(finalScore2{e}==0) = 0;
        MT = [MT;tmp(:)];
        P = [P;[tmp3(:) tmp2(:) e*ones(size(tmp(:)))]]; 
    end
end
plot3(P(:,1),P(:,2),P(:,3),'.')
title([num2str(mean(abs(MT))) '--' num2str(corr(P(:,1),P(:,2)))]);
kidx = P(:,1) ~= 0 & P(:,2) ~= 0;
mean(abs(MT(kidx)))
hold on
plot3(P(kidx,1),P(kidx,2),P(kidx,3),'ro')
%% apply BEST
%[score] = detectEmergence(FileList{5},emergenceNet,BESTEmergence,151,zU,zE,11,[],[]);
%[score1,score2,raw1,raw2] = detectEmergence(FileList{1},emergenceNet3,emergenceDynamic,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
[finalScore] = detectEmergence(FileList{2},'',emergenceDynamic,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
%[score1,score2,raw1,raw2] = detectEmergence(FileList{1},'',emergenceDynamic,151,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
%%
func = @(X)detectEmergence(X,emergenceNet3,emergenceDynamic,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
pf = partialFunction(func,'cornPopperNNapp');
pf.publish();
%%
layers = [imageInputLayer([size(rZ,1) size(rZ,2) 3]);
          convolution2dLayer([12 12],20);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer([5 5],20);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
layers = [imageInputLayer([size(rZ,1) size(rZ,2) 3]);
          convolution2dLayer([21 21],7);
          reluLayer();
          maxPooling2dLayer(8,'Stride',2);
          %convolution2dLayer([5 5],20);
          %reluLayer();
          %maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
layers = [imageInputLayer([size(rZ,1) size(rZ,2) 3]);
          convolution2dLayer([11 11],9);
          reluLayer();
          maxPooling2dLayer([8 8],'Stride',6);
          convolution2dLayer([3 3],4);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
%selIdx = randperm(numel(TYPE));
%selIdx = 1:249:numel(TYPE);
%subX = imgSIM(:,:,:,1:4:end);
%subY = TYPE(1:4:end);
options = trainingOptions('sgdm','MaxEpochs',5,'InitialLearnRate',0.0001,'ExecutionEnvironment','parallel');
emergenceNet3 = trainNetwork(single(rZ),categorical(TYPE),layers,options);
%emergenceNet4 = trainNetwork(subX,categorical(subY),layers,options);
%%
TESTN = predict(emergenceNet4,imgSIM(:,:,:,1:3*249));
%%
[score] = detectEmergence(FileList{3},emergenceNet,151);
%% regress 
layers = [ ...
    imageInputLayer([size(zCfinal,1) size(zCfinal,2) 1])
    convolution2dLayer([2 15],10)
    reluLayer
    maxPooling2dLayer([2 2],'Stride',1);
    %convolution2dLayer([10 1],5)
    %reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','MaxEpochs',700,'InitialLearnRate',0.00001,'ExecutionEnvironment','parallel');
nVALUE = VALUE;
nVALUE(find(isinf(VALUE))) = 0;
UnVALUE = mean(nVALUE);
nVALUE = nVALUE - UnVALUE;
RtrainedNet = trainNetwork(single(zCfinal),nVALUE,layers,options);

%% TEST BEST
close all
[YBEST PROB] = classify(BESTEmergence,single(BEST));
plot(YBEST)
%% test prediction for NON regression
TEST = predict(emergenceNet,single(rZ));
%% test regression
close all
RTEST = predict(RtrainedNet,single(zCfinal)) + UnVALUE;
plot(VALUE,RTEST,'.')
%% reshape
nsz = 101;
rZ = zeros([nsz nsz 3 numel(Z)]);
for e = 1:numel(Z)
    rZ(:,:,:,e) = Z{e};
    e
end
%% check a test
Tnum = 2;
[Tpth] = fileparts(FileList{Tnum}{1});
fidx = strfind(Tpth,filesep);
% make labels
LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
HL = {'d' 'p'};
LABELS = {};
for e1 = 1:numel(HL)
    for e2 = 1:numel(LL)
        for e3 = 1:numel(NL)
            LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
        end
    end
end
LABELS = lower(LABELS);
TESTValues = [];
for l = 1:numel(LABELS)
    Tkey{l} = lower([Tpth((fidx(end)+1):end) '-' LABELS{l}]);
    sc = handData_t2.get(Tkey{l});
    TESTValues = [TESTValues ; sc];
end
TESTValues = reshape(TESTValues,[250 168]);
Tscore = [];
for t = 1:size(TESTValues,2)
    fidx = find(TESTValues(:,t));
    if ~isempty(fidx)
        Tscore(t) = fidx(1);
    else
        Tscore(t) = 0;
    end
end
%%
DELTA = abs(double(score2) - Tscore);
[~,sidx] = sort(DELTA);
%%

%%
Rscore = [];
mag = 5;
for t = 1:size(score,1)
    fidx = find(squeeze(score(t,2,:) > .5));
    if ~isempty(fidx)
        Rscore(t) = mag*fidx(1);
    else
        Rscore(t) = -1;
    end
end
close all
DELTA = abs(Tscore - Rscore);
[J,sidx] = sort(DELTA);
figure;
plot(Tscore,Rscore,'.')
hold on
plot(Tscore(sidx(end)),Rscore(sidx(end)),'ro');
figure;
plot(squeeze(score(sidx(end),2,:)))
fidx = find(Rscore == 10);
figure;
plot(squeeze(score(fidx(1),2,:)))

%%
        
        
        SZ = size(STACK);
        STACK = reshape(STACK,[prod(SZ(1:2)) SZ(3:5)]);
        STACK = permute(STACK,[1 3 2 4]);
        
        
        
        close all
        for cb = 1:size(STACK,5)
            for t = 1:size(STACK,4)
                imshow(STACK(:,:,:,t,cb),[]);
                drawnow
            end
        end
        
        

%%
close all
for e = 1:size(CheckerBoard,1)
    [cameraParameters,imagePoints,worldPoints] = getCameraPara(CheckerBoard{e,1}{1});
    for s = 1:100
        tic
       
        [rI{e}(:,:,s)] = getRectifiedImage(CheckerBoard{e,1}{s},cameraParameters,imagePoints,worldPoints);
        toc
    end
    
    %{
    SZ = size(subI{e});
    tform = fitgeotrans(imagePoints,worldPoints,'affine');
    
    
    
    cI = imwarp(subI{e},tform);
    
    
    
    E = edge(rgb2gray(subI{e}),'Canny');
    [H,T,R] = hough(E','theta',linspace(-15,15,200));
    P  = houghpeaks(H,40,'NHoodSize',[101 11]);
    lines = houghlines(E',T,R,P,'FillGap',600,'MinLength',400);
    
    [H2,T2,R2] = hough(E,'theta',linspace(-20,20,200));
    P2  = houghpeaks(H2,40,'NHoodSize',[101 11]);
    lines2 = houghlines(E,T2,R2,P2,'FillGap',600,'MinLength',400);
    
    close all
    figure, imshow(subI{e}), hold on
    max_len = 0;
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       xy = flipdim(xy,2);
       
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    
    max_len = 0;
    for k = 1:length(lines2)
       xy = [lines2(k).point1; lines2(k).point2];
       
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines2(k).point1 - lines2(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    
    
    plot(imagePoints(:,1),imagePoints(:,2),'k*');
    
    %{
    imshow(subI{e},[]);
    drawnow
    %}
    drawnow
    %}
end

%% find frist images
for e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    nm = str2num(nm);
    if nm == 1
        kp(e) = true;
    else
        kp(e) = false;
    end
end
FileList = FileList(kp);
%%
dlP = '/home/nate/Downloads/JEFF/';
mkdir(dlP);
cnt = 1;
for e = 30:numel(FileList)
    CMD = ['iget ' FileList{e} ' ' dlP num2str(cnt) '.tif'];
    cnt = cnt + 1;
    system(CMD)
    CMD
end
%%
FilePath = 'W:\';
FilePath = '/home/nate/Downloads/JEFF/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1)
%%
I = {};
for e = 1:numel(FileList)
    try
        I{end+1} = imread(FileList{e});
        imshow(I{end},[]);
        drawnow
    catch
    end
    
end
%%
un = CheckerBoard.checkerBoard;
for e = 1:size(un,1)
    CB{e} = un(e,:);
end
CheckerBoard.checkerBoard = CB'
%%
options = trainingOptions('sgdm','MaxEpochs', 10,'ExecutionEnvironment','auto');%,'InitialLearnRate', 1e-6);
%layers = data.layers;
layers = [imageInputLayer([75 75 3])
          convolution2dLayer([15,15],10)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          convolution2dLayer([7,7],3)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer()
          classificationLayer()];


checkerDetector = trainFasterRCNNObjectDetector(CheckerBoard, layers, options,'SmallestImageDimension',400,'NumStrongestRegions',10);


%%
img = imread(CheckerBoard{1,1}{1});
[bbox, score, label] = detect(checkerDetector, img,'SelectStrongest',false);

detectedImg = insertShape(img, 'Rectangle', bbox);
figure
imshow(detectedImg)


%% get values for a
