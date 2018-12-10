%% super train
FilePath = '/mnt/tetra/nate/bitPile/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
load(['/mnt/tetra/nate/bitPile/prepData.mat'],'FE','FU');
Y = {};
X = {};
XCNN = {};

parfor e = 1:numel(FileList)
    
    loadResult = load(FileList{e})
    [FC,toFeatureFORCNN] = prepareForNetworks(loadResult.sig,FE,FU);
    
    
    [pth,nm,ext] = fileparts(FileList{e});
    fidx = strfind(pth,filesep);
    type = pth((fidx(end)+1):end);
    
    if strcmp(type,'1')
        %Y = [Y;1*ones(size(FC,1),1)];
        Y{e} = [1*ones(size(FC,1),1)];
    else
        %Y = [Y;0*ones(size(FC,1),1)];
        Y{e} = [0*ones(size(FC,1),1)];
    end
    
    fsz = size(toFeatureFORCNN);
    toFeatureFORCNN = reshape(toFeatureFORCNN,[prod(fsz(1:3)) fsz(4)]);
    
    X{e} = FC;
    XCNN{e} = toFeatureFORCNN';
    %X = [X;FC];
    %XCNN = cat(4,XCNN,toFeatureFORCNN);
    
    e
end
%%
nX = cell2mat(X');
nY = cell2mat(Y');
cnnX = cell2mat(XCNN');
%%
newnewfda = fitcdiscr(newnewC,nY,'DiscrimType','quadratic');
newnewNET = patternnet([11]);
newnewNET = train(newnewNET,newnewC',full(ind2vec((nY+1)')),'useParallel','yes');
%%
szX = size(cnnX);
%
layers = [imageInputLayer(szX(1:3));
          convolution2dLayer(5,13);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
      options = trainingOptions('sgdm','Plots','training-progress','MaxEpochs',20,'InitialLearnRate',0.0001,'ExecutionEnvironment','parallel');

      newCNN = trainNetwork(cnnX,categorical(nY),layers,options);
%% gather all sorghumData
dataPath = ['/iplant/home/leakey_cyverse/sorghumData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        sorFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%% try whole fft
oI = imread(sorFileList{2});
%% try the repeated smooth and subtract
I = oI;
for r = 1:5
    background = imfilter(I,fspecial('gaussian',[101 101],51),'replicate');
    I = I - background;
    imshow(background,[])
    drawnow
end
figure;
imshow(I,[]);
%% reset I
I = oI;
%% perform ops
F = fft2(I - mean(I(:)));
amp = abs(F);
ang = angle(F);
[sortAmp,sidx] = sort(amp(:),'descend');
%% generate fft image
newAmp = zeros(size(amp));
numF = 2;
newAmp(sidx(1:numF)) = sortAmp(1:numF);
newF = newAmp.*exp(i*ang);
look = ifft2(newF);
close all
figure;imshow(I,[]);
figure;imshow(look,[]);
%% gather all maizeData
dataPath = ['/iplant/home/leakey_cyverse/maizeData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        maizeFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%% gather all setariaData
dataPath = ['/iplant/home/leakey_cyverse/setariaData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        setFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%%
toH = [];
for e = 1:20
    tmpH = imread(sorFileList{e});
    toH = [toH;tmpH(:)];
    imshow(tmpH,[]);
    e
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scan for training data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = ['/iplant/home/phytomorphuser/workIT%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
fPatchFileList = {};
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.mat')
        fPatchFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
fPatchFileList
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull training set local
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
massDownload('/iplant/home/phytomorphuser/workIT/', '.mat','/mnt/tetra/nate/workIT/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scan local for mat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cMatFilePath = '/mnt/tetra/nate/workIT/';
fPatchFileList = {};
FileExt = {'mat'};
fPatchFileList = gdig(cMatFilePath,fPatchFileList,FileExt,1);
%% build out basis vector - from CONDOR
[basisU,basisE] = buildGradeBasis(fPatchFileList,100);
%% test

    %[data] = myRemoteLoader(gradeFileList{1},'T');
%data = prepareData(fPatchFileList{1},basisU,basisE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build out human click map(s) - from CONDOR
%%%%%%%%%%%%%%%%%%%%%%%%
close all
%%%%%%%%%%%%%%%%%%%%%%%%
% generate a temp local location
%%%%%%%%%%%%%%%%%%%%%%%%
tmpPath = ['/mnt/tetra/nate/tmp/' tempname filesep];
mkdir(tmpPath);
humanClickPath = '/mnt/tetra/nate/sorStomataClick/';
mkdir(humanClickPath);
toClick = 100;
for e = 1:toClick
    try
        tmpM = zeros(size(rI));
        [pth,nm,ext] = fileparts(fPatchFileList{e});
        outFile = [humanClickPath nm '.mat'];
        fileName = fPatchFileList{e};
        if ~exist(outFile)
            fileList{1} = fileName;
            %fileList = xfer_get({fileName},tmpPath,0,0);
            load(fileList{1},'T','rI','customData');

            tmpM(customData{2}) = 1;
            R = regionprops(logical(tmpM),'boundingBox');
            subI = imcrop(rI,R(1).BoundingBox);
            uix = [];
            [uix(:,2),uix(:,1),~] = impixel(subI,[]);

            imshow(subI,[]);
            hold on
            plot(uix(:,2),uix(:,1),'r*')
            hold off
            drawnow
            save(outFile,'T','rI','customData','uix');
        else
            e
        end
    catch ME
        ME
    end
end
%% build out X and Y for training - from CONDOR
close all
humanClickPath = '/mnt/tetra/nate/sorStomataClick/';
trainFileList = {};
FileExt = {'mat'};
trainFileList = gdig(humanClickPath,trainFileList,FileExt,1);
dilateValue = 9;
[labeledTrainingPackage] = generateMaskPackage(trainFileList,100,dilateValue,basisU,basisE);
%% train labeled package - CONDOR workflow
[AI_layer] = trainAIlayer(labeledTrainingPackage,[1 1]);
%% applyAI layer - CONDOR
[probMap] = applyAIlayer(fPatchFileList{1},AI_layer,basisU,basisE);
%% optimize AI output
[opti_para] = optimizeAIoutput(AI_layer,labeledTrainingPackage);
%% test opti
[probMap] = applyAIlayer(fPatchFileList{4},AI_layer,basisU,basisE);
[grade,ret1,ret2,BO] = opti(probMap,'',opti_para,labeledTrainingPackage.oI(1:108,1:108,4));
%% test full package
%{
%% execute extractor CONDOR -- #1
oPath = '/mnt/tetra/nate/bitPile/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read example image and get expected image size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpI = imread(sorFileList{ss});
tmpI = imhistmatch(tmpI,toH);
%tmpI = tmpI(end/4:3*end/4,end/4:3*end/4);
%tmpI = STORE(:,:,1);
imageSZ = size(tmpI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate operational grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pointSet,pointSetSize] = genIXgrid(imageSZ,[1 1],[40 40]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sample domain-1 disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opDomain = {};
domainSize = {};
R = [0 40];
N = [41 250];
[n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
X = n1.*cos(n2+pi);
Y = n1.*sin(n2+pi);
opDomain{1} = [X(:) Y(:)];
domainSize{1} = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sample domain-2 disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[X,Y] = ndgrid(-30:30,-30:30);
%opDomain{2} = [X(:) Y(:)];
%domainSize{2} = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init feature extraction and start sample/extract phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NF = 35;
N = [41 250];
eFunc = @(X)extractFFTfeatures(X,NF,N);
[T,rI] = ixExtractor(tmpI,tmpI,R,pointSet,opDomain,domainSize,eFunc,false);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate training data on CONDOR for sorghum data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate operational grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borderSize = [40 40];
[pointSet,pointSetSize] = genIXgrid(imageSZ,[1 1],borderSize);
innerSize = imageSZ - 2*borderSize;
numMini = [4 4];
miniSZ = ceil(innerSize.*numMini.^-1);
miniLabel = [floor((pointSet(:,1)-borderSize(1)-1).*miniSZ(1).^-1) floor((pointSet(:,2)-borderSize(2)-1).*miniSZ(2).^-1)];
[UQ,i1,i2] = unique(miniLabel,'rows');
subPointSet = {};
for u = 1:size(UQ,1)
    subPointSet{u} = pointSet(i2==u,:);
    globalIDX{u} = find(i2==u);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sample domain-1 disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opDomain = {};
domainSize = {};
R = [0 40];
N = [41 250];
[n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
X = n1.*cos(n2+pi);
Y = n1.*sin(n2+pi);
opDomain{1} = [X(:) Y(:)];
domainSize{1} = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sample domain-2 disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[X,Y] = ndgrid(-30:30,-30:30);
%opDomain{2} = [X(:) Y(:)];
%domainSize{2} = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package input(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trimValues = R;
pointSetData = pointSet;
domainData{1} = opDomain;
domainData{2} = domainSize;
domainData{3} = trimValues;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issue ticket(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TsorFileList = issueBulkTicket(sorFileList);
rPath = '/iplant/home/nmiller/workIT/';
rPath = '/iplant/home/phytomorphuser/workITOUT/';
[rPath iticket] = issueTicket(rPath(1:end-1),100*numel(TsorFileList),'write');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issue ticket(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = './outputFeatures/';
fibratedExtrationLayer = cFlow('extractionLayer');
fibratedExtrationLayer.setMCRversion('v930');
for e = 1:300
    for sD = 1%:numel(subPointSet)
        
        [pth,nm,ext] = fileparts(sorFileList{e});
        customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(sD) '}'];
        customData{2} = globalIDX{sD};
        
        % for local run
        % extractionLayer(sorFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
        fibratedExtrationLayer(TsorFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
        %[forSHOWEX,fMASKEX] = applyAllLayers(sorFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,pointSet,customData,AI_layer,basisU,basisE,opti_para,oPath,rPath);
        %cnt = cnt + 1;
        sD
    end
    e
end
fibratedExtrationLayer.submitDag(auth,150,150);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issue ticket(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = './outputFeatures/';
func = cFlow('applyAllLayers');
func.setMCRversion('v930');
for e = 1
    for sD = 1%:numel(subPointSet)
        
        [pth,nm,ext] = fileparts(sorFileList{e});
        customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(sD) '}'];
        customData{2} = globalIDX{sD};
        
        % for local run
        func(TsorFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,pointSet,customData,AI_layer,basisU,basisE,opti_para,oPath,rPath);
       
        sD
    end
    e
end
func.submitDag(auth,150,150);
    
%%
e = 4;
rPath = '';
oPath = '';
[forSHOWEX,fMASKEX] = applyAllLayers(sorFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,pointSet,customData,AI_layer,basisU,basisE,opti_para,oPath,rPath);

%% select project - try
C = [];
FL = {};
newnewC = [];
for e = 1:size(cnnX,3)
    
    tmp = squeeze(cnnX(:,:,e,:));
   
    
    
    if e ~= 1
        tmp = cumsum(tmp,2);
        %tmp = diff(tmp,1,4);
        %tmp = diff(tmp,1,5);
    end
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:2)) sz(3)])';
    
    selectIDX = find(nY==1);
    selectIDX0 = find(nY==0);
    
    [FUSEL{e},FESEL{e},FL{e}] = PCA_FIT_FULLws(tmp(selectIDX,:),15);
    
    FCSEL{e} = PCA_REPROJ(tmp,FESEL{e},FUSEL{e});
    SIM = PCA_BKPROJ(FCSEL{e},FESEL{e},FUSEL{e});
    
    DEL = SIM - tmp;
    DEL = sum(DEL.*DEL,2).^.5;
    
    %FC{e} = reshape(FC{e},[sz(1:3) 10]);
    %FC{e} = permute(FC{e},[1 2 4 3]);
    newnewC = cat(2,newnewC,[FCSEL{e} DEL]);
end
%%
close all
pileOutput = '/mnt/tetra/nate/bitPile/';
mkdir([pileOutput '0' filesep]);
mkdir([pileOutput '1' filesep]);
for ss = 1:50%numel(sorFileList)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read example image and get expected image size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpI = imread(sorFileList{ss});
    tmpI = imhistmatch(tmpI,toH);
    %tmpI = tmpI(end/4:3*end/4,end/4:3*end/4);
    %tmpI = STORE(:,:,1);
    imageSZ = size(tmpI);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate operational grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pointSet,pointSetSize] = genIXgrid(imageSZ,[1 1],[40 40]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate sample domain-1 disk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opDomain = {};
    domainSize = {};
    R = [0 40];
    N = [41 250];
    [n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
    X = n1.*cos(n2+pi);
    Y = n1.*sin(n2+pi);
    opDomain{1} = [X(:) Y(:)];
    domainSize{1} = size(X);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate sample domain-2 disk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[X,Y] = ndgrid(-30:30,-30:30);
    %opDomain{2} = [X(:) Y(:)];
    %domainSize{2} = size(X);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init feature extraction and start sample/extract phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NF = 35;
    N = [41 250];
    eFunc = @(X)extractFFTfeatures(X,NF,N);
    [T,rI] = ixExtractor(tmpI,tmpI,R,pointSet,opDomain,domainSize,eFunc,false);
    
    
    
    
    dataIndex = [4 5] - 3;
    featureImageSize = size(rI);
    opData1 = thawTensor(T,dataIndex(1));
    opData2 = thawTensor(T,dataIndex(2));
    toFeature = cat(4,opData1,cos(opData2),sin(opData2));
  
    toFeature = permute(toFeature,[2 1 4 3]);
    toFeatureFORCNN = toFeature;
    toFeature(:,:,2,:) = cumsum(toFeature(:,:,2,:),2);
    toFeature(:,:,3,:) = cumsum(toFeature(:,:,3,:),2);
    %toFeatureFORCNN = toFeature;
    
    
    
    toFeature = permute(toFeature,[1 2 4 3]);
    FC = [];
    for d = 1:size(toFeature,4)
        tmp = toFeature(:,:,:,d);
        sz = size(tmp);
        tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);
        tmpFC = PCA_REPROJ(tmp',FESEL{d},FUSEL{d});
        
        
        
        
        
        SIM = PCA_BKPROJ(tmpFC,FESEL{e},FUSEL{e});
        DEL = SIM - tmp;
        DEL = sum(DEL.*DEL,2).^.5;
    
        
        
        FC = cat(2,FC,[tmpFC DEL]);
    end
    % z-score normalize
    %FC = bsxfun(@minus,FC,mu);
    %FC = bsxfun(@times,FC,sigma.^-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start prediction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convolution network predict
    CNNY = newCNN.predict(toFeatureFORCNN);
    CNNY = double(CNNY(:,2));
    %%%%%%%%%%%%%%%%%%%%
    % predict with linear discriminate analysis
    %%%%%%%%%%%%%%%%%%%%
    [yPRE,yPROB] = predict(newfda,FC);
    %yPROB_G = sum(yPROB(:,clusterN(gSEL,1)+1:end),2);
    yPROB_G = sum(yPROB(:,2:end),2);
    %%%%%%%%%%%%%%%%%%%%
    % predict with neural network
    %%%%%%%%%%%%%%%%%%%%
    yPROB2 = newNET(FC')';
    %yPROB = sum(yPROB2(:,clusterN(gSEL,1)+1:end),2);
    yPROB = sum(yPROB2(:,2:end),2);
    %%%%%%%%%%%%%%%%%%%%
    % stack prediction
    %%%%%%%%%%%%%%%%%%%%
    PTOT = reshape([yPROB_G yPROB CNNY],[featureImageSize(1),featureImageSize(2) 3]);
    % filter prediction
    %PTOT = imfilter(PTOT,fspecial('gaussian',[21 21],5),'replicate');
    %%%%%%%%%%%%%%%%%%%%
    % execute optimized operations for stomata isolation 
    %%%%%%%%%%%%%%%%%%%%
    [qqq,a,outX,fMASK,pcount] = opti(PTOT,[],g5,rI,init,delta);

    %{
        pout = nSTORE(:,:,:,1); 
[qqq,a,outX,fMASK,TC] = opti(pout,[],g5,rI,init,delta);
    %}
    %{
    fidx0 = find(fMASK==0);
    fidx1 = find(fMASK==1);
    
    fidx0 = fidx0(randperm(numel(fidx0)));
    fidx0 = fidx0(1:numel(fidx1));
    
    
    sig = T(:,fidx0);
    save([pileOutput '0' filesep num2str(ss) '.mat'],'sig');
    sig = T(:,fidx1);
    save([pileOutput '1' filesep num2str(ss) '.mat'],'sig');
    %}
    imshow(outX,[]);
    title(num2str(ss));
    drawnow
    nSTORE(:,:,:,ss) = outX;
    mSTORE(:,:,:,ss) = PTOT;
    iSTORE(:,:,ss) = rI;
end

%% run local over sub data
clear retVec1 retVec2
R = [0 40];
N = [41 250];
NF = 35;
mag = 1;
disp = false;
func = [];
NUM = 100;
tempy = imread(sorFileList{1});
tempy = tempy(end/4:3*end/4,end/4:3*end/4);
%tempy = tempy(end/2 - 40:end/2 + 40,end/2 - 40:end/2 + 40);
[TretVec1,TretVec2,tOI] = gogo_bugEye(tempy,tempy,R,N,NF,mag,[],1,disp,func);
retVec1 = zeros([size(TretVec1) NUM]);
%retVec2 = zeros([size(TretVec2) NUM],'single');
oI = zeros([size(tOI) NUM]);
clear TrevVec1 TretVec2
rot = linspace(-pi,pi,13);
cnt = 1;

cnt = 1;
for e = 1:NUM
    tic
    rawT = imread(sorFileList{e});
    rawT = imhistmatch(rawT,toH);
    %{
    for rep = 1:3
        BK = imfilter(rawT,fspecial('gaussian',[171 171],81),'replicate');
        %BK = imfilter(rawT,fspecial('disk',[3ave1]),'replicate');
        rawT = rawT - BK;
    end
    rawT = bindVec(rawT);
    %}
    
    
    %rawT(1,:) = [];
    %rawT(:,1) = [];

    %for r = 1:numel(rot)
        %rawI = imrotate(rawT,rot(r)*180/pi,'bicubic','crop');
        %rawI = rawT((end-1)/2 - 40:(end-1)/2 + 40,(end-1)/2 - 40:(end-1)/2 + 40);
        %rawI = imrotate(rawI,rot(r)*180/pi,'bicubic','crop');
        %rawI = rawI((end-1)/2 - 40:(end-1)/2 + 40,(end-1)/2 - 40:(end-1)/2 + 40);
        %imshow(rawI,[]);
        %drawnow
        rawI = rawT(end/4:3*end/4,end/4:3*end/4);
        %rawI = imrotate(rawT,rot(r)*180/pi,'bicubic','crop');
        %toOp = rawI;
        toOp = rawI;
        STORE(:,:,cnt) = rawI;
        %[retVec1(:,:,:,:,:,cnt),retVec2(:,:,:,:,cnt),oI(:,:,cnt)] = gogo_bugEye(rawI,toOp,R,N,NF,mag,[],1,disp,func);
        [retVec1(:,:,:,:,:,cnt),~,oI(:,:,cnt)] = gogo_bugEye(rawI,toOp,R,N,NF,mag,[],1,disp,func);
        
        cnt = cnt + 1;
    %end
    toc
    e
end
%%
retVec2(isnan(retVec2(:)) | isinf(retVec2(:))) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASTERTHIS
%%
toOp = permute(retVec1,[1 2 6 3 4 5]);
clear retVec1
toOp = cat(6,toOp(:,:,:,:,:,1),cos(toOp(:,:,:,:,:,2)),sin(toOp(:,:,:,:,:,2)));
toOp = permute(toOp,[1 2 3 5 4 6]);
%%

%%
C = [];
FL = {};
FC = {};
for e = 1:size(toOp,6)
    tmp = toOpF(:,:,:,:,:,e);
   
    if e ~= 1
        %tmp = diff(tmp,1,4);
        %tmp = diff(tmp,1,5);
    end
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:3)) prod(sz(4:5))]);
    [FU{e},FE{e},FL{e}] = PCA_FIT_FULLws(tmp,10);
    FC{e} = PCA_REPROJ(tmp,FE{e},FU{e});
    FC{e} = reshape(FC{e},[sz(1:3) 10]);
    FC{e} = permute(FC{e},[1 2 4 3]);
    C = cat(3,C,FC{e});
end
%%
tryClickAgain(oI,C);
%%
cp = {};
for e = 1:size(oI,3)
    [cp{e}(:,2),cp{e}(:,1),~] = impixel(oI(:,:,e),[]);
    e
end
%% for each cp - EXTRACT FROM USERS CLICKS
close all
uC = [];
uL = [];
for e = 1:numel(cp)
    tmpM = zeros(size(oI,1),size(oI,2));
    tmpF = C(:,:,:,e);
    szC = size(tmpF);
    tmpF = reshape(tmpF,[prod(szC(1:2)) szC(3)]);
    tmpL = zeros(size(tmpF,1),1);
    
    for p = 1:size(cp{e},1)
        IDX = sub2ind([size(oI,1) size(oI,2)],cp{e}(:,1),cp{e}(:,2));
        tmpM(IDX) = 1;
    end
    
    tmpM = imdilate(tmpM,strel('disk',11,0));
    IDX = find(tmpM==1);
    out = flattenMaskOverlay(bindVec(oI(:,:,e)),logical(tmpM));
    imshow(out,[]);
    drawnow
    uC = [uC;tmpF];
    uL = [uL;tmpM(:)];
    sMSK(:,:,e) = tmpM;
end
%%
NET = feedforwardnet(13);
NET = patternnet([7]);
NET = train(NET,nX',full(ind2vec((nL)')),'useParallel','yes');
%% sample equal parts from each group
UQ = unique(YT);
%UQ = unique(uL);
NM = [];
for u = 1:numel(UQ)
    NM(u) = sum(YT==UQ(u));
    %NM(u) = sum(uL==UQ(u));
end
mm = min(NM);
mag = [1 1 1 1 2 2];
nX = [];
nL = [];
for e = 1:numel(UQ)
    fidx = find(YT==UQ(e));
    sub = X(fidx,:);
    fidx = randperm(size(sub,1));
    nX = [nX;sub(fidx(1:mag(e)*mm),:)];
    nL = [nL;e*ones(mag(e)*mm,1)];
end
%%
ClassTreeEns = fitensemble(X,uL);
%%
SVM = fitcsvm(X(1:100:end,1:5),uL(1:100:end));
%%
[IDX, Z] = rankfeatures(nX(1:10:end,:)',nL(1:10:end)','Criterion','entropy','NumberOfIndices',10);
%%
Mdl = fitcdiscr(nX(:,IDX),nL,'DiscrimType','quadratic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GOODNESS
%% zero layer
%[X,mu,sigma] = zscore(uC);
X = uC;
UQL = unique(uL);
clusterN = [[1 1];[1 2];[1 3];[1 4];[2 1];[2 2];[2 3];[2 4]];
clusterN = [1 1];
Y_LABELS = [];
XSAM = {};
YSAM = {};
for cl = 1:size(clusterN,1)
    
    YT = uL;
    for u = 1:numel(UQL)
        if u == 1
            sam = 10;
        else
            sam = 1;
        end
        clear GMModel_sub;
        F = find(uL == UQL(u));
        
        options = statset('Display','iter','MaxIter',300);
        GMModel_sub = fitgmdist(X(F(1:sam:end),:),clusterN(cl,u),'Options',options,'Replicates',3);
        [subIDX] = GMModel_sub.cluster(X(F,:));
        if u > 1
            subIDX = subIDX + clusterN(cl,u-1);
        end
        YT(F) = subIDX;
    end
    
    Y_LABELS = [Y_LABELS YT];
    
end

for cl = 1:size(clusterN,1)
    YT = Y_LABELS(:,cl);
    UQ = unique(YT);
    
    NM = [];
    for u = 1:numel(UQ)
        NM(u) = sum(YT==UQ(u));
    end
    
    mm = min(NM);
    mag = [1 1 1 1 1 1];
    nX = [];
    nL = [];
    for e = 1:numel(UQ)
        fidx = find(YT==UQ(e));
        sub = X(fidx,:);
        fidx = randperm(size(sub,1));
        
        nX = [nX;sub(fidx(1:mag(e)*mm),:)];
        %nX = [nX;sub];
        
        nL = [nL;e*ones(mag(e)*mm,1)];
        %nL = [nL;e*ones(size(sub,1),1)];
    end
    
    
    XSAM{cl} = nX;
    YSAM{cl} = nL;
end
%% zero 2
clear marginStats L
for cl = 1:size(clusterN,1)
    FLDA{cl} = fitcdiscr(XSAM{cl},YSAM{cl},'DiscrimType','quadratic','CrossVal','on');
    L{cl} = kfoldLoss(FLDA{cl});
    [label,score,cost] = kfoldPredict(FLDA{cl});
    m = kfoldMargin(FLDA{cl});
    marginStats{cl} = table(min(m),mean(m),max(m),'VariableNames',{'Min','Mean','Max'});
    marginStats{cl} 
end
%%
gSEL = 1;
fda = fitcdiscr(XSAM{gSEL},YSAM{gSEL},'DiscrimType','quadratic');
NET = patternnet([11]);
NET = train(NET,XSAM{gSEL}',full(ind2vec((YSAM{gSEL})')),'useParallel','yes');
%%
szC = size(toOpF);
szC(3) = 50;
cnnX = reshape(toOpF(:,:,1:50,:,:,:),[prod(szC(1:3)) szC(4:6)]);
cnnX = permute(cnnX,[2 3 4 1]);
%
layers = [imageInputLayer([szC(4:6)]);
          convolution2dLayer(5,13);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
      options = trainingOptions('sgdm','Plots','training-progress','MaxEpochs',20,'InitialLearnRate',0.0001,'ExecutionEnvironment','parallel');

      CNN = trainNetwork(cnnX(:,:,:,1:5:end),categorical(uL(1:5:end)),layers,options);
   %%
   
szC = size(toOpF);
cnnX = reshape(toOpF(:,:,:,:,:,:),[prod(szC(1:2)) prod(szC(3)) szC(4:6)]);

%%
COLLECT = true;

if COLLECT
    SS = [];
    cTH = [];
    
    CNTF = true;
    disp = false;
else
    
    CNTF = false;
    if numel(cTH) == 5
        disp = true;
    end
end
close all
h1 = figure;
h2 = figure;



for e = 1:50%size(C,4)%numel(cp)%
    tmpF = C(:,:,:,e);
   
    szC = size(tmpF);
    
    TT = squeeze(cnnX(:,e,:,:,:));
    TT = permute(TT,[2 3 4 1]);
    CNNY = CNN.predict(TT);
    CNNY = double(CNNY(:,2));
   
    
    
    
    tmpF = reshape(tmpF,[prod(szC(1:2)) szC(3)]);
    %tmpF = bsxfun(@minus,tmpF,mu);
    %tmpF = bsxfun(@times,tmpF,sigma.^-1);
    
   
    
    [yPRE,yPROB] = predict(fda,tmpF);
    yPROB_G = sum(yPROB(:,clusterN(gSEL,1)+1:end),2);
    
    
    yPROB2 = NET(tmpF')';
    yPROB = sum(yPROB2(:,clusterN(gSEL,1)+1:end),2);
    
    
    yPROB1 = reshape(yPROB_G,[size(oI,1),size(oI,2)]);
    yPROB2 = reshape(yPROB,[size(oI,1),size(oI,2)]);
    
    yPROBM = prod([yPROB_G yPROB CNNY],2);
    yPROB = mean([yPROB_G yPROB CNNY],2);
    sPROB = std([yPROB_G yPROB CNNY],1,2);
    
    PSTACK(:,:,:,e) = reshape([yPROB_G yPROB CNNY],[size(oI,1),size(oI,2) 3]);
    
    CNNY = reshape(CNNY,[size(oI,1),size(oI,2)]);
    sPROB = reshape(sPROB,[size(oI,1),size(oI,2)]);
    yPROB = reshape(yPROB,[size(oI,1),size(oI,2)]);
    yPROBM = reshape(yPROBM,[size(oI,1),size(oI,2)]);
    
    %yPROB = yPROB + 1*sPROB;
    
    MASKM = yPROBM > .5^3;
    MASKS = yPROB > .5;
    
    
    
    
    STOREYY(:,:,e) = yPROBM;
    imshow(oI(:,:,e),[]);
    marker = imfilter(yPROBM,fspecial('gaussian',[21 21],3),'replicate');
    %marker = interp2(marker,2);
    hold on
    DS = [];
    cnt = 1;
    l = 1;
    imgSZ = size(oI);
    for LEVEL = linspace(.009,.2,30)
        
        mm = marker - LEVEL;
        RECON = imreconstruct(mm,marker);
        BO = (mm - RECON) > - LEVEL;
        
        
        R = regionprops(BO,yPROB,'Centroid','PixelIdxList','Area','MeanIntensity','Eccentricity','MajorAxisLength');
        for r = 1:numel(R)
            DS(cnt).Area = R(r).Area;
            DS(cnt).PixelIdxList = R(r).PixelIdxList;
            DS(cnt).MeanIntensity = R(r).MeanIntensity;
            DS(cnt).Eccentricity = 1 - R(r).Eccentricity;
            DS(cnt).MajorAxisLength = R(r).MajorAxisLength;
            delta = (min(abs(R(r).Centroid(1) - imgSZ(2))).^2 + min(abs(R(r).Centroid(2) - imgSZ(1))).^2).^.5;
            DS(cnt).DistanceToEdge = delta;
            DS(cnt).level = l;
            cnt = cnt + 1;
        end
        
        
        
        
        dB = bwboundaries(BO);
        hold on
        for b = 1:numel(dB)
            plot(dB{b}(:,2),dB{b}(:,1),'r')
            hold on
        end
        l = l + 1;
    end
    
    LEVELS = max([DS.level]);
    
    for l = 1:LEVELS
        lidx = find([DS.level] == l);
        for a = 1:numel(lidx)
            childIDX = DS(lidx(a)).PixelIdxList;
            childLEVEL = DS(lidx(a)).level;
            for p = 1:numel(DS)
                if p ~= lidx(a) && (DS(p).level) == childLEVEL + 1
                    inter = intersect(DS(p).PixelIdxList,childIDX);
                    if ~isempty(inter)
                        DS(lidx(a)).parent = p;
                        break
                    end
                end
            end
            fprintf(['Done searching: ' num2str(a) ':' num2str(numel(lidx)) '@level' num2str(l) '\n']);
        end
    end
    
    topList = find([DS.level] == 1);
    FL = [[DS.Area]' [DS.Eccentricity]' [DS.MeanIntensity]'];
    filter = find((FL(:,2) > .25) & (FL(:,3) > .3) & FL(:,1) > 6^2*pi);
    myMask1 = zeros(size(MASK));
    for o = 1:numel(filter)
        myMask1(DS(filter(o)).PixelIdxList) = 1;
    end
    
    RS = regionprops(logical(myMask1),'MajorAxisLength');
    singleLength = mean([RS.MajorAxisLength]);
    
    
    
    myMask2 = zeros(size(myMask1));
    filter = find(round([DS.MajorAxisLength]' * singleLength^-1)==2 & (FL(:,3) > .6));
    for o = 1:numel(filter)
        myMask2(DS(filter(o)).PixelIdxList) = 1;
    end
    
    for r = 1:4
        orginalBorder{r} = myMask1(:,1);
        myMask1(:,1) = 1;
        myMask1 = imrotate(myMask1,90);
    end
    holeFilled = imfill(myMask1,'holes');
    holeFilled = myMask1 == 0 & holeFilled == 1;
    holeFilledBK = bwlarge(holeFilled);
    holeFilled = holeFilledBK == 0 & holeFilled == 1;
    myMask1 = myMask1 + holeFilled;
    for r = 1:4
        myMask1(:,1) = orginalBorder{r};
        myMask1 = imrotate(myMask1,90);
    end
    
    for r = 1:4
        orginalBorder{r} = myMask2(:,1);
        myMask2(:,1) = 1;
        myMask2 = imrotate(myMask2,90);
    end
    holeFilled = imfill(myMask2,'holes');
    holeFilled = myMask2 == 0 & holeFilled == 1;
    holeFilledBK = bwlarge(holeFilled);
    holeFilled = holeFilledBK == 0 & holeFilled == 1;
    myMask2 = myMask2 + holeFilled;
    for r = 1:4
        myMask2(:,1) = orginalBorder{r};
        myMask2 = imrotate(myMask2,90);
    end
    
    % upgrade single to doubles
    RS1 = regionprops(logical(myMask1),'PixelIdxList');
    RS2 = regionprops(logical(myMask2),'PixelIdxList');
    myMask1 = zeros(size(myMask1));
    myMask2 = zeros(size(myMask2));
    
    toRM = [];
    for r = 1:numel(RS1)
        for s = 1:numel(RS2)
            if ~isempty(intersect(RS2(s).PixelIdxList,RS1(r).PixelIdxList))
                toRM = [toRM r];
            end
        end
    end
    RS1(toRM) = [];
    for r = 1:numel(RS1)
        myMask1(RS1(r).PixelIdxList) = 1;
    end
    for r = 1:numel(RS2)
        myMask2(RS2(r).PixelIdxList) = 1;
    end
    
    out = flattenMaskOverlay(oI(:,:,e),logical(myMask1),.5,'r');
    out = flattenMaskOverlay(out,logical(myMask2),.5,'g');
    if e <= 50
        out = flattenMaskOverlay(out,logical(sMSK(:,:,e)),.5,'b');
    end
    imshow(out,[],'Border','tight');
    drawnow
    RRO(:,:,:,e) = out;
    DSSS{e} = DS;
    
    %{
    %%%%%%%%%%%%%%%%
    thresholdValues = [6^2*pi .25 .3 15 6^2*pi .25 .3];
    thresholdValues = [xn3 15 xn2];
    thresholdValues = xn5;
    func = @(X)generateStomataMask(DSSS,sMSK,X);
    ops = optimset('Display','iter');
    
    xn5 = fminsearch(func,thresholdValues,ops);
    [grade,myMask1,myMask2,outx3,agHOPE2] = generateStomataMask(DSSS,sMSK,g5,oI);
    [grade,myMask1,myMask2,out3,ag2] = generateStomataMask(DSSS,sMSK,thresholdValues,oI);
    
    
    for e = 1:size(out2,4)
        %TOP = [repmat(oI(:,:,e),[1 1 3]) out2(:,:,:,e) outx(:,:,:,e)];
        BOTTOM = [repmat(oI(:,:,e),[1 1 3]) outX(:,:,:,e)];
        %BOTTOM = [repmat(oI(:,:,e),[1 1 3]) out2(:,:,:,e) outx2(:,:,:,e)];
        imshow([BOTTOM],[]);
        hold on
        plot(cp{e}(:,2),cp{e}(:,1),'r*');
        title(num2str(e));
        
        %[newX,newY,V] = impixel();
        %plot(newX,newY,'g*');
        %cp{e} = [cp{e};[newY newX]];
        drawnow
        waitforbuttonpress
        hold off
    end
    
    alpha = [10 100 100 100 2 10 2 1 10];
    
    delta = [0.1 0.1 0.1 0.1 10 0.1 1- 1 1]
    
    para = [.1 1*3/9 1*3/9 7*3/9 80 .3 30 3 3];
    
    init = [.1 1*3/9 1*3/9 7*3/9 100 .3 30 3 3];
    delta = [.05 .2 .2 .2 10 0.15 10 1 1]*.05^-1;
    %para = [g1 50 3];
    %para(5) = 100;
    para = para.*alpha;
    func = @(X)opti(PSTACK,sMSK,X,oI,init,delta);
    ops = optimset('Display','iter');
    g5new = fminsearch(func,g5,ops);
    [qqq,a,outX] = opti(PSTACK,sMSK,g5,oI,init,delta);
    
    
    %{
    
    for c = 1:numel(topList)
        M(c).Area = [];
        M(c).m0 = [];
        curIDX = topList(c);
        while ~isempty(DS(curIDX).parent)
            M(c).Area = [M(c).Area DS(curIDX).Area];
            M(c).m0 = [M(c).m0 DS(curIDX).Area*DS(curIDX).Eccentricity*DS(curIDX).MeanIntensity];
            curIDX = DS(curIDX).parent;
        end
        plot(M(c).m0)
        hold on
        drawnow
    end
    R = regionprops(logical(myMask1),'PixelIdxList','Area');
    cnt = count([R.Area]);
    R(cnt ~= 1) = [];
    myMask1 = zeros(size(myMask1));
    for r = 1:numel(R)
        myMask1(R(r).PixelIdxList) = 1;
    end
    %}
    
    
    
    
    
    
    
    %imshow(BO,[]);
    %drawnow
    break
    
    MASKS = bwareaopen(MASKS,round(5^2*pi));
    
    RECON = imreconstruct(MASKS,BO);
    MASK = MASKS & BO; 
    MASK = imclose(MASK,strel('disk',4,0));
    
    %break
   
    
    
    
    
    
    
    if CNTF
        R = regionprops(MASK,yPROB,'PixelValues','PixelList','Area','PixelIdxList','BoundingBox','Eccentricity');
        cidx0 = count([R.Area]);
        R(cidx0 ~= 1 ) = [];
        MASK = zeros(size(MASK));
        for c = 1:numel(R)
            MASK(R(c).PixelIdxList) = 1;
        end
        MASK = logical(MASK);
    end
    
    
    
    
    
    
    out = flattenMaskOverlay(oI(:,:,e),MASK);
    
    if COLLECT
        
        msk = zeros(size(msk));
        for c = 1:numel(R)
            msk(R(c).PixelIdxList) = 1;
            tmpB = imcrop(yPROB,R(c).BoundingBox);
            tmpB = imresize(tmpB,[18 18]);
            SS = cat(3,SS,tmpB);
        end
        
        
    else
        
        filt = mean(SS,3);
        filt = filt / norm(filt(:));
        func = @(X)X(:)'*filt(:)/norm(X(:));
        B = nlfilter(yPROB,size(filt),func);
        
        
        
        %{
        marker = imfilter(B,fspecial('gaussian',[21 21],7),'replicate');
        RECON = imreconstruct(marker,B);
        BO = ~(RECON == B);
        RO = regionprops(BO,yPROB,'MeanIntensity','PixelIdxList','MaxIntensity');
        out = flattenMaskOverlay(oI(:,:,e),BO);
        imshow(out,[]);
        drawnow
        %}
        
        BO = imdilate(B,strel('disk',11,0));
        BO = BO == B;
        BO = imdilate(BO,strel('disk',11,0));
        RO = regionprops(BO,yPROB,'MeanIntensity','PixelIdxList','MaxIntensity');
        
        
        
        if e <= 50
            here = 1
            if numel(cTH) ~= 50
                what = 1
                FV = linspace(.3,.7,10);
                TP = [];
                FP = [];

                for f = 1:numel(FV)
                    MM = zeros(size(yPROB));
                    TOT(f) = 0;


                    
                    
                    tmpM = sMSK(:,:,e);

                    for b = 1:numel(RO)
                        if RO(b).MeanIntensity > FV(f)
                            MM(RO(b).PixelIdxList) = 1;
                            TOT(f) = TOT(f) + 1;
                            if tmpM(RO(b).PixelIdxList) == 0
                                FP(f,b) = 1;
                            end
                        else
                            FP(f,b) = 0;
                        end
                    end
                    
                   
                    


                    for b = 1:size(cp{e},1)
                        TP(f,b) = MM(cp{e}(b,1),cp{e}(b,2));
                    end



                end

                FN = sum(TP'==0,1);
                FP = sum(FP'==1,1);
                TC = (FP + FN);
                [v,TH(e)] = min(TC);
                cTH(e) = mean(FV(find(TC == v)));
            end
            
            
            
            if numel(cTH) == 50
                
                MM = zeros(size(yPROB));
                for b = 1:numel(RO)
                    if RO(b).MeanIntensity > mean(cTH)
                        MM(RO(b).PixelIdxList) = 1;
                        TOT(f) = TOT(f) + 1;
                    end
                end
                MASK = logical(MM);
                
                 
                MM = zeros(size(yPROB));
                for b = 1:numel(RO)
                    if RO(b).MeanIntensity > (mean(cTH) + std(cTH))
                        MM(RO(b).PixelIdxList) = 1;
                        TOT(f) = TOT(f) + 1;
                    end
                end
                gMASK = logical(MM);
                
                 
                MM = zeros(size(yPROB));
                for b = 1:numel(RO)
                    if RO(b).MeanIntensity > (mean(cTH) - std(cTH))
                        MM(RO(b).PixelIdxList) = 1;
                        TOT(f) = TOT(f) + 1;
                    end
                end
                lMASK = logical(MM);
                
                out = flattenMaskOverlay(oI(:,:,e),lMASK,.3,'r');
                out = flattenMaskOverlay(out,MASK,.3,'r');
                out = flattenMaskOverlay(out,gMASK,.3,'r');
                %out = flattenMaskOverlay(oI(:,:,e),MASK,.3,'r');
            end
            
        end
        
        
        
        
    end
    
    
    
    
    
    figure(h1);
    imshow([out repmat(yPROB,[1 1 3]) repmat(yPROB1,[1 1 3]) repmat(yPROB2,[1 1 3]) repmat(CNNY,[1 1 3])],[]);
    title(num2str(numel(cTH)))
    
    
    if COLLECT
        figure(h2);
        plot(cTH);
    else
        if e == 1
            close(h2);
        end
    end
    drawnow
    if disp
        waitforbuttonpress
    end
    %}
end
%% zero view

close all
TYPES = [];
COLLECT = false;


if COLLECT
    SS = [];
end


for e = 1:size(C,4)%numel(cp)%
    tmpF = C(:,:,:,e);
   
    szC = size(tmpF);
    tmpF = reshape(tmpF,[prod(szC(1:2)) szC(3)]);
    tmpF = bsxfun(@minus,tmpF,mu);
    tmpF = bsxfun(@times,tmpF,sigma.^-1);
    
    
    
    [yPRE,yPROB] = predict(fda,tmpF);
    yPROB_G = yPROB(:,clusterN(gSEL,1)+1:end);
    
    
    yPROB2 = NET(tmpF')';
    yPROB = sum(yPROB2(:,clusterN(gSEL,1)+1:end),2);
    
    yPROB = mean([yPROB_G yPROB],2);
    
    
    
    
    
    %yPROB = NET(tmpF')';
    %{
    WHAT = NET(tmpF');
    WHAT = WHAT(2,:);
    %}
    
    %RGB = yPROB(:,1:3);
    
    %yPROB = sum(yPROB(:,1:cluster1),2);
    
    %yPROB = reshape(yPROB,[size(oI,1),size(oI,2)]);
    %yPROB = .5*(yPROB + yPROB2);
    
    
    yPROB = reshape(yPROB,[size(oI,1),size(oI,2)]); 
    
    
    
    
    RGB = reshape(RGB,[size(oI,1),size(oI,2) 3]);
    %yPROB = imfilter(yPROB,fspecial('gaussian',[31 31],4));
    
    
    
    if ~COLLECT
        filt = mean(SS,3);
        filt = filt / norm(filt(:));
        %FF = xcorr2(yPROB,filt);
        %yPROB2 = imfilter(yPROB,filt,'replicate');
        func = @(X)X(:)'*filt(:)/norm(X(:));
        B = nlfilter(yPROB,size(filt),func);


        BO = imdilate(B,strel('disk',11,0));
        BO = BO == B;
        BO = imdilate(BO,strel('disk',9,0));

        RO = regionprops(BO,yPROB,'MeanIntensity','PixelIdxList');
    end
    
    
    if COLLECT
        if e <= 50
            FV = linspace(.3,.7,10);
            TP = [];
            for f = 1:numel(FV)
                MM = zeros(size(yPROB));
                TOT(f) = 0;

                for b = 1:numel(RO)
                    if RO(b).MeanIntensity > FV(f)
                        MM(RO(b).PixelIdxList) = 1;
                        TOT(f) = TOT(f) + 1;
                    end
                end


                for b = 1:size(cp{e},1)
                    TP(f,b) = MM(cp{e}(b,1),cp{e}(b,2));
                end



            end

            FN = sum(TP'==0,1);
            FP = TOT - size(cp{e},1);
            TC = (FP + FN);
            [~,TH(e)] = min(TC);
        end
    end

    
    
    MM = zeros(size(yPROB));
    for b = 1:numel(RO)
        if RO(b).MeanIntensity > FV(round(mean(TH)))
            MM(RO(b).PixelIdxList) = 1;
            TOT(f) = TOT(f) + 1;
        end
    end
    
    
    
    
    msk = BO;
    if ~COLLECT
        msk = MM;
    else
        msk = yPROB > .5;
    end
    pp = find(msk == 1);
    msk = bwareaopen(msk,50);
    
    
    %{
    [yPRE2,yPROB2] = predict(MODEL2,tmpF(pp,:));
    yPROB2 = sum(yPROB2(:,1:cluster12),2);
    yPROB2M = zeros(size(yPROB));
    yPROB2M(pp) = yPROB2(:);
    %}
    %{
    fp = msk & ~logical(sMSK(:,:,e));
    fn = ~msk & logical(sMSK(:,:,e));
    tp = msk & logical(sMSK(:,:,e));
    tn = ~msk & ~logical(sMSK(:,:,e));
    %}
    
    
    prod(size(msk))
    newL = [0*ones(sum(fp(:)),1);1*ones(sum(tp(:)),1);2*ones(sum(fn(:)),1);3*ones(sum(tn(:)),1)];
    numel(newL)
    TYPES = [TYPES;newL];
    
    R = regionprops(msk,yPROB,'PixelValues','PixelList','Area','PixelIdxList','BoundingBox');
    
    cidx0 = count([R.Area]);
    R(cidx0 ~= 1) = [];
    
    if COLLECT
        msk = zeros(size(msk));
        for c = 1:numel(R)
            msk(R(c).PixelIdxList) = 1;
            tmpB = imcrop(yPROB,R(c).BoundingBox);
            tmpB = imresize(tmpB,[18 18]);
            SS = cat(3,SS,tmpB);
        end
    end
   
    
    
    
    
    R2 = regionprops(msk,yPROB2M,'PixelValues','PixelList','MeanIntensity');
    dB = bwboundaries(msk);
    imshow([repmat(bindVec(oI(:,:,e)),[1 1 3]) repmat(yPROB,[1 1 3]) RGB repmat(msk,[1 1 3])],[0 1]);
    hold on
    for r = 1:numel(dB)
        if true%R2(r).MeanIntensity > .45
            plot(dB{r}(:,2),dB{r}(:,1),'g');
        else
            plot(dB{r}(:,2),dB{r}(:,1),'r');
        end
        hold on
    end
    for r = 1:numel(R)
        [~,l] = max(R(r).PixelValues);
        plot(R(r).PixelList(l,1),R(r).PixelList(l,2),'r.');
    end
    
    if e < 50
        plot(cp{e}(:,2),cp{e}(:,1),'g*');
    end
    

    
    
    hold off
    title(num2str(e))
    drawnow
    pause(.004)
    if e > 0
        waitforbuttonpress
    end
end
%% first layer
[X,mu,sigma] = zscore(uC);
%X = uC;


cluster1 = 1;
clear GMModel_sub;
F1 = find(uL == 1);
options = statset('Display','iter','MaxIter',300);
GMModel_sub = fitgmdist(X(F1,:),cluster1,'Options',options,'Replicates',3);
[subIDX1] = GMModel_sub.cluster(X(F1,:));

cluster2 = 3;
F0 = find(uL==0);
sam = 10;
GMModel_sub = fitgmdist(X(F0(1:sam:end),:),cluster2,'Options',options,'Replicates',3);
[subIDX2] = GMModel_sub.cluster(X(F0,:));

YT = uL;
YT(F1) = subIDX1;
YT(F0) = subIDX2 + cluster1;
MODEL = fitcnb(X,YT);
%% second layer
cluster12 = 2;
clear GMModel_sub;
F1 = find(TYPES == 1);
options = statset('Display','iter','MaxIter',300);
GMModel_sub = fitgmdist(X(F1,:),cluster12,'Options',options,'Replicates',3);
[subIDX1] = GMModel_sub.cluster(X(F1,:));

cluster22 = 2;
F0 = find(TYPES == 0);
GMModel_sub = fitgmdist(X(F0,:),cluster22,'Options',options,'Replicates',3);
[subIDX2] = GMModel_sub.cluster(X(F0,:));

YT2(F1) = subIDX1;
YT2(F0) = subIDX2 + cluster12;
YT2 = [subIDX1;subIDX2+cluster12];
MODEL2 = fitcnb(X([F1;F0],:),YT2);
%%
close all
[lambda] = myLDA(X(:,:),uL);
cLL = X(:,:)*lambda;
ksdensity(cLL(uL==0));
hold on
ksdensity(cLL(uL==1));
hold on
%%
close all
[lambda] = mynLDA(X(:,6:end),uL,1,4);
cLL = X(:,6:end)*lambda;

for e = 1:size(cLL,2)
    figure;
    ksdensity(cLL(uL==0,e));
    hold on
    ksdensity(cLL(uL==1,e));
    hold on
end
%% ksdensity
for comp = 1:15
    close all
    [f,xi] = ksdensity(X(F1,comp));
    plot(xi,f,'r')
    hold all
    [f,xi] = ksdensity(X(F0,comp));
    plot(xi,f,'k')
    hold all
    CLL = {'r--' 'r--' 'r--' 'b--' 'b--' 'b--' 'b--' 'b--' 'b--' 'b--' 'b--' 'b--'};
    for e = 1:numel(unique(YT))
        [f,xi] = ksdensity(X(YT==e,comp));
        plot(xi,f,CLL{e})
    end
    waitforbuttonpress
end
%%
close all
TYPES = [];
%SS = [];
for e = 1:size(C,4)%numel(cp)%
    tmpF = C(:,:,:,e);
   
    szC = size(tmpF);
    tmpF = reshape(tmpF,[prod(szC(1:2)) szC(3)]);
    tmpF = bsxfun(@minus,tmpF,mu);
    tmpF = bsxfun(@times,tmpF,sigma.^-1);
    [yPRE,yPROB] = predict(Mdl,tmpF);
    yPROB_G = yPROB(:,1);
    
    
    yPROB2 = NET(tmpF')';
    yPROB = sum(yPROB2(:,1),2);
    
    yPROB = mean([yPROB_G yPROB],2);
    
    
    
    
    
    %yPROB = NET(tmpF')';
    %{
    WHAT = NET(tmpF');
    WHAT = WHAT(2,:);
    %}
    
    %RGB = yPROB(:,1:3);
    
    %yPROB = sum(yPROB(:,1:cluster1),2);
    
    %yPROB = reshape(yPROB,[size(oI,1),size(oI,2)]);
    %yPROB = .5*(yPROB + yPROB2);
    
    
    yPROB = reshape(yPROB,[size(oI,1),size(oI,2)]); 
    
    
    
    
    RGB = reshape(RGB,[size(oI,1),size(oI,2) 3]);
    %yPROB = imfilter(yPROB,fspecial('gaussian',[31 31],4));
    
    filt = mean(SS,3);
    filt = filt / norm(filt(:));
    %FF = xcorr2(yPROB,filt);
    %yPROB2 = imfilter(yPROB,filt,'replicate');
    func = @(X)X(:)'*filt(:)/norm(X(:));
    B = nlfilter(yPROB,size(filt),func);
   
    
    BO = imdilate(B,strel('disk',11,0));
    BO = BO == B;
    BO = imdilate(BO,strel('disk',9,0));
    
    RO = regionprops(BO,yPROB,'MeanIntensity','PixelIdxList');
   
    
    if e <= 50
        FV = linspace(.3,.7,10);
        TP = [];
        for f = 1:numel(FV)
            MM = zeros(size(yPROB));
            TOT(f) = 0;

            for b = 1:numel(RO)
                if RO(b).MeanIntensity > FV(f)
                    MM(RO(b).PixelIdxList) = 1;
                    TOT(f) = TOT(f) + 1;
                end
            end


            for b = 1:size(cp{e},1)
                TP(f,b) = MM(cp{e}(b,1),cp{e}(b,2));
            end



        end

        FN = sum(TP'==0,1);
        FP = TOT - size(cp{e},1);
        TC = (FP + FN);
        [~,TH(e)] = min(TC);
    end

    
    
    MM = zeros(size(yPROB));
    for b = 1:numel(RO)
        if RO(b).MeanIntensity > FV(round(mean(TH)))
            MM(RO(b).PixelIdxList) = 1;
            TOT(f) = TOT(f) + 1;
        end
    end
    
    
    
    msk = yPROB > .5 & BO;
    msk = BO;
    msk = MM;
    pp = find(msk == 1);
    msk = bwareaopen(msk,50);
    
    
    
    [yPRE2,yPROB2] = predict(MODEL2,tmpF(pp,:));
    yPROB2 = sum(yPROB2(:,1:cluster12),2);
    yPROB2M = zeros(size(yPROB));
    yPROB2M(pp) = yPROB2(:);
    %{
    fp = msk & ~logical(sMSK(:,:,e));
    fn = ~msk & logical(sMSK(:,:,e));
    tp = msk & logical(sMSK(:,:,e));
    tn = ~msk & ~logical(sMSK(:,:,e));
    %}
    
    
    prod(size(msk))
    newL = [0*ones(sum(fp(:)),1);1*ones(sum(tp(:)),1);2*ones(sum(fn(:)),1);3*ones(sum(tn(:)),1)];
    numel(newL)
    TYPES = [TYPES;newL];
    
    R = regionprops(msk,yPROB,'PixelValues','PixelList','Area','PixelIdxList','BoundingBox');
    
    cidx0 = count([R.Area]);
    R(cidx0 ~= 1) = [];
    
    %{
    msk = zeros(size(msk));
    for c = 1:numel(R)
        msk(R(c).PixelIdxList) = 1;
        tmpB = imcrop(yPROB,R(c).BoundingBox);
        tmpB = imresize(tmpB,[18 18]);
        SS = cat(3,SS,tmpB);
    end
    %}
   
    
    
    
    
    R2 = regionprops(msk,yPROB2M,'PixelValues','PixelList','MeanIntensity');
    dB = bwboundaries(msk);
    imshow([repmat(bindVec(oI(:,:,e)),[1 1 3]) repmat(yPROB,[1 1 3]) RGB repmat(msk,[1 1 3])],[0 1]);
    hold on
    for r = 1:numel(dB)
        if true%R2(r).MeanIntensity > .45
            plot(dB{r}(:,2),dB{r}(:,1),'g');
        else
            plot(dB{r}(:,2),dB{r}(:,1),'r');
        end
        hold on
    end
    for r = 1:numel(R)
        [~,l] = max(R(r).PixelValues);
        plot(R(r).PixelList(l,1),R(r).PixelList(l,2),'r.');
    end
    
    if e < 50
        plot(cp{e}(:,2),cp{e}(:,1),'g*');
    end
    

    
    
    hold off
    title(num2str(e))
    drawnow
    pause(.004)
    if e > 0
        %waitforbuttonpress
    end
end
%% 
figure;
for e = 1:size(C,4)
    for k = 1:size(C,3)
        imshow(C(:,:,k,e),[]);
        drawnow
        pause(.1)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try strand per K for retVec1 - slice 1
numPC = 15;
rSZ = size(retVec1);
strandBundleB = squeeze(retVec1(:,:,:,:,1,:));
strandBundleB = permute(strandBundleB,[1 2 4 5 3]);
rSZW = size(strandBundleB);
SIMY1 = [];
for kk = 1:rSZW(3)
    strandBundleBT = reshape(squeeze(strandBundleB(:,:,kk,:,:)),[prod(rSZW([1 2 4])) rSZW(end)]);
    [Uk1{kk},Ek1{kk},Lk1{kk}] = PCA_FIT_FULLws(strandBundleBT,numPC);
    Ck1{kk} = PCA_REPROJ(strandBundleBT,Ek1{kk},Uk1{kk});
    SIMk1{kk} = PCA_BKPROJ(Ck1{kk},Ek1{kk},Uk1{kk});
    SIMY1 = [SIMY1;SIMk1{kk}];
    kk
end
%% try strand per K for retVec1 - slice 2
close all
numPC = 15;
rSZ = size(retVec1);
strandBundleB = squeeze(retVec1(:,:,:,:,2,:));
strandBundleB = permute(strandBundleB,[1 2 4 5 3]);
rSZW = size(strandBundleB);
SIMY2 = [];
for kk = 1:rSZW(3)
    strandBundleBT = reshape(squeeze(strandBundleB(:,:,kk,:,:)),[prod(rSZW([1 2 4])) rSZW(end)]);
    [Uk2{kk},Ek2{kk},Lk2{kk}] = PCA_FIT_FULLws(strandBundleBT,numPC);
    Ck2{kk} = PCA_REPROJ(strandBundleBT,Ek2{kk},Uk2{kk});
    SIMk2{kk} = PCA_BKPROJ(Ck2{kk},Ek2{kk},Uk2{kk});
    SIMY2 = [SIMY2;SIMk2{kk}];
    kk
    hope = fft(strandBundleBT,[],2);
    plot(mean(abs(hope),1))
    hold all;
    drawnow
    H(:,kk,:) = abs(hope(:,2:10));
end
%% try strand per K for retVec2
close all
numPC = 15;
rSZ = size(retVec2);
strandBundleB = real(squeeze(retVec2(:,:,:,:,:)));
strandBundleB = permute(strandBundleB,[1 2 5 4 3]);
rSZW = size(strandBundleB);
SIMY2 = [];
for kk = 1:rSZW(3)
    strandBundleBT = reshape(squeeze(strandBundleB(:,:,:,kk,:)),[prod(rSZW([1 2 4])) rSZW(end)]);
    [Uk2{kk},Ek2{kk},Lk2{kk}] = PCA_FIT_FULLws(strandBundleBT,numPC);
    Ck2{kk} = PCA_REPROJ(strandBundleBT,Ek2{kk},Uk2{kk});
    SIMk2{kk} = PCA_BKPROJ(Ck2{kk},Ek2{kk},Uk2{kk});
    SIMY2 = [SIMY2;SIMk2{kk}];
    kk
    hope = fft(strandBundleBT,[],2);
    plot(mean(abs(hope),1))
    hold all;
    drawnow
    H(:,kk,:) = abs(hope(:,2:10));
end
%%
%%%%%%%%%%%%%%%%%%%%%%
clear SIMk2 SIMk1
SIMY1 = reshape(SIMY1,[rSZW(1:2) rSZW(4) rSZW(3) rSZW(5)]);
SIMY1 = permute(SIMY1,[1 2 5 4 3]);
rSZ = size(SIMY1);
SIMY1 = reshape(SIMY1,[rSZ(1:4) 1 rSZ(5)]);

SIMY2 = reshape(SIMY2,[rSZW(1:2) rSZW(4) rSZW(3) rSZW(5)]);
SIMY2 = permute(SIMY2,[1 2 5 4 3]);
rSZ = size(SIMY2);
SIMY2 = reshape(SIMY2,[rSZ(1:4) 1 rSZ(5)]);

SIMYT = cat(5,SIMY1,SIMY2);
%% watch sims per strand
strandBundleB = squeeze(retVec1(:,:,:,:,2,:));
strandBundleB = permute(strandBundleB,[1 2 4 5 3]);
close all
for kk = 1:rSZW(3)
    strandBundleBT = reshape(squeeze(strandBundleB(:,:,kk,:,:)),[prod(rSZW([1 2 4])) rSZW(end)]);
    for s = 1:100:size(strandBundleBT,2)
        plot(strandBundleBT(s,:)');
        hold on
        plot(SIMk2{kk}(s,:)')
        hold off
        drawnow
    end
end

%% cluster based on k = k0
k0 = [2 3];
CLUS = 6;
sam = 100;
tmpKM = [];
fI = [];
for k = 1:numel(k0)
    tmpKM = [tmpKM squeeze(reshape(Ck2{k0(k)},[prod(rSZW([1:2 4])) 1 numPC]))];
    W = squeeze(reshape(Ck2{k0(k)},[(rSZW([1 2 4])) 1 numPC]));
    fI = cat(3,fI,permute(W,[1 2 4 3]));
end
%%
CLUS = 7;
tmpKM = reshape(H(:,:,:),[size(H,1) size(H,2)*size(H,3)]);
options = statset('Display','iter','MaxIter',300);
GMModel = fitgmdist(tmpKM(1:sam:end,:),CLUS,'Options',options,'RegularizationValue',0.001);
kCluster = GMModel.cluster(tmpKM);
kCluster = reshape(kCluster,[rSZW([1:2 4])]);
% view k0 cluster
close all
for e = 1:size(kCluster,3)
    RGB = label2rgb(kCluster(:,:,e));
    imshow([double(RGB)/255 repmat(bindVec(oI(:,:,e)),[1 1 3])],[]);
    drawnow
    waitforbuttonpress
end
%%
close all
tryClickAgain(oI,fI,retVec1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% try strand analysis
%close all
numPC = 20;

rSZ = size(retVec1);
strandBundleB = squeeze(retVec1(:,:,:,:,1,:));
strandBundleB = permute(strandBundleB,[1 2 4 5 3]);
rSZW = size(strandBundleB);
strandBundleB = reshape(strandBundleB,[prod(rSZW(1:4)) rSZW(end)]);
[Ua1,Ea1,La1] = PCA_FIT_FULLws(strandBundleB,numPC);
Ca1 = PCA_REPROJ(strandBundleB,Ea1,Ua1);
Sa1 = PCA_BKPROJ(Ca1,Ea1,Ua1);

rSZ = size(retVec1);
strandBundleB = squeeze(retVec1(:,:,:,:,2,:));
initPHI = strandBundleB(:,:,1,:,:);
strandBundleB = bsxfun(@minus,strandBundleB,initPHI);
strandBundleB = permute(strandBundleB,[1 2 4 5 3]);
rSZW = size(strandBundleB);
strandBundleB = reshape(strandBundleB,[prod(rSZW(1:4)) rSZW(end)]);
[Ua2,Ea2,La2] = PCA_FIT_FULLws(strandBundleB,numPC);
Ca2 = PCA_REPROJ(strandBundleB,Ea2,Ua2);
Sa2 = PCA_BKPROJ(Ca2,Ea2,Ua2);
%% stack simulated strands
Sa1 = reshape(Sa1,rSZW);
Sa1 = permute(Sa1,[1 2 5 3 4]);
rSZ = size(Sa1);
Sa1 = reshape(Sa1,[rSZ(1:4) 1 rSZ(5)]);

Sa2 = reshape(Sa2,rSZW);
Sa2 = permute(Sa2,[1 2 5 3 4]);
rSZ = size(Sa2);
Sa2 = reshape(Sa2,[rSZ(1:4) 1 rSZ(5)]);

SaT = cat(5,Sa1,Sa2);
initPHI = reshape(initPHI,[size(initPHI,1) size(initPHI,2) size(initPHI,3) size(initPHI,4) 1 size(initPHI,5)]);
%% STRAND
MS = [];
CC = [];
fSZ = size(SaT);
TOCOMPRESS = 200;
for level = 1:2
    toOp = squeeze(SaT(:,:,:,:,level,:));
    tmpSZ = size(toOp);
    masterStack = reshape(toOp,[tmpSZ(1:2) prod(tmpSZ(3:4)) tmpSZ(end)]);
    masterStack = permute(masterStack,[1 2 4 3]);
    mSZ = size(masterStack);
    masterStack = reshape(masterStack,[prod(mSZ(1:3)) mSZ(4)]);
    [Um{level},Em{level},Lm{level}] = PCA_FIT_FULLws(masterStack,TOCOMPRESS);
    Cm = PCA_REPROJ(masterStack,Em{level},Um{level});
    
    masterSIM = PCA_BKPROJ(Cm,Em{level},Um{level});
    masterSIM = reshape(masterSIM,mSZ);
    masterSIM = ipermute(masterSIM,[1 2 4 3]);
    masterSIM = reshape(masterSIM,tmpSZ);
    masterSIM = reshape(masterSIM,[tmpSZ(1:4) 1 tmpSZ(5)]);
    MS = cat(5,MS,masterSIM);
    mSZ(end) = TOCOMPRESS;
    Cm = reshape(Cm,mSZ);
    Cm = ipermute(Cm,[1 2 4 3]);
    
    Cm = reshape(Cm,[mSZ(1:2) mSZ(4) 1 mSZ(3)]);
    CC = cat(4,CC,Cm);
    humm = 1
end
%CC = permute(CC,[1 2 4 5 3]);
%%
close all
for e = 3:size(CC,5)
    for s = 1:size(CC,4)
        ZZ = zeros(size(oI,1),size(oI,2));
        WH = ZZ;
        for l = 1:size(CC,3)
            tempy = CC(:,:,l,s,e);
            
            tempy = tempy - mean(tempy(:));
            zz = fft2(tempy);
            ampy = abs(fftshift(zz));
            ed = edge(tempy);
            ampy = abs(zz);
            [sA,iA] = sort(ampy(:),'descend');
            ampy = imfilter(ampy,fspecial('gaussian',[21 21],3),'replicate');
            rA = ampy;
            ampy = bindVec(ampy);
           
            msk = ampy > 1*graythresh(ampy);
            ampy = ifftshift(msk .* rA);
            ampy = zeros(size(ampy));
            ampy(iA(1:30)) = sA(1:30);
            iF = ifft2(ampy.*angle(-i*zz));
            
            
            [dx] = corner(tempy);
            IDX = sub2ind(size(ZZ),dx(:,2),dx(:,1));
            %ZZ(IDX) = ZZ(IDX) + 1;
            ZZ = ZZ + iF;
            WH = WH + ampy;
            out = flattenMaskOverlay(bindVec(oI(:,:,e)),ed);
            %imshow([out repmat(bindVec(ZZ),[1 1 3]) repmat(bindVec(iF),[1 1 3]) repmat(bindVec(tempy),[1 1 3])  repmat(bindVec(abs(fftshift(zz))),[1 1 3])],[],'Border','tight');
            imshow([out repmat(bindVec(tempy),[1 1 3])  repmat(bindVec(abs(fftshift(zz))),[1 1 3])],[],'Border','tight');
            
            hold on
            plot(dx(:,1),dx(:,2),'g.')
            title(num2str(l))
            drawnow
            %waitforbuttonpress
            pause(.02)
            hold off
        end
        waitforbuttonpress
    end
end
%%
h1= figure;
h2= figure;
for sample = 1:size(Cm,4)

    for toView = 1:10
        figure(h1)
        imshow(bindVec(Cm(:,:,toView,sample)),[])
        figure(h2)
        imshow(oI(:,:,sample),[]);
        drawnow
        waitforbuttonpress
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
sam = 1000;

kS = squeeze(strandBundleB(1:sam:end,:,3,:));
toC = repmat(1:size(kS,2),[size(kS,1) 1]);
kS = cat(3,toC,kS);
figure;
hold on
initPhi = kS(:,1,3);
kS(:,:,3) = bsxfun(@minus,kS(:,:,3),initPhi);
tmp = kS;
tmp = permute(tmp,[1 3 2]);
tmp3 = squeeze(tmp(:,3,:));
[Sa1 Ca Ua Ea La ERRa LAMa] = PCA_FIT_FULL(tmp3,10);

%%
szP = size(tmp);
for e = 2:3
    kS(:,:,e) = bindVec(kS(:,:,e));
end
for e = 1:300%size(kS,1)
    plot3(squeeze(kS(e,:,1)),squeeze(kS(e,:,2)),squeeze(kS(e,:,3)),'r');
    drawnow
end
%% whole wave analysis for patch strip pull down
close all
for e = 1:10:100
    close all
    I = imread(sorFileList{e});

    for rep = 1:3
        BK = imfilter(I,fspecial('gaussian',[171 171],81),'replicate');
        %BK = imfilter(rawT,fspecial('disk',[31]),'replicate');
        I = I - BK;
    end

    F = fft2(I - mean(I(:)));
    amp = abs(F);
    ang = angle(F);
    [sortAmp,sidx] = sort(amp(:),'descend');
    Z = zeros(size(amp));

    
    [k1,k2] = ind2sub(size(I),sidx);
    mA = atan2(size(I,2)*(k2-1).^-1,size(I,1)*(k1-1).^-1);
    
    %{
    close all
    plot(mA(1:200),sortAmp(1:200),'.')
    waitforbuttonpress
    %}
    
    options = statset('Display','iter','MaxIter',300);
    GMModel = fitgmdist(mA(1:100),3,'Options',options,'RegularizationValue',0.001);
    cidx = GMModel.cluster(mA(1:100));
    
    DIR = GMModel.mu;
    %{
    for u = 1:3
        
        fidx = find(cidx==u);
        
        [ampU,ms] = sort(sortAmp(fidx),'descend');
        
        Z(sidx(fidx(ms(1:2)))) = ampU(1:2);%sortAmp(1:2);
    end
    %}
    %{
    
    hi = histogram(mA,1000);
    [~,midx] = max(hi.BinCounts);
    hidx = mA > hi.BinEdges(midx) & mA < hi.BinEdges(midx+1);
    DIR = mean(mA(hidx));
    %}
    
   
    
    %DIR = mean(mA);
    %DIR = median(mA);
    %{
    HIX = linspace(hi.BinLimits(1),hi.BinLimits(2),10);

    DIR = HIX(midx);
    %}
    VEC = [cos(DIR) sin(DIR)];

   
    simI = ifft2(Z.*exp(-1i*ang));
    out = flattenMaskOverlay(bindVec(I),bindVec(simI) > .5,.4);
     imshow(out,[])
     hold on
     %figure;
     %imshow(simI,[]);
    %hold on
    %imshow(simI,[]);
    %figure;
    for e = 1:size(VEC,1)
        quiver(size(I,2)/2,size(I,1)/2,VEC(e,1),VEC(e,2),100,'LineWidth',5);
        title(num2str(e))
    end
    waitforbuttonpress
end
%% wave sim
%close all
PZ = [R(2) R(2)];
[n1,n2] = ndgrid(-PZ(1):PZ(1),-PZ(2):PZ(2));
rho = (n1.^2 + n2.^2).^.5;
theta = atan2(n1,n2);
PT = [104 104];
PT = [50 70];
waveT = squeeze(retVec1(PT(1),PT(2),:,:,:,1));
%waveT = squeeze(SaT(PT(1),PT(2),:,:,:,1));
%waveT = squeeze(MS(PT(1),PT(2),:,:,:,1));
%waveT = squeeze(SIMYT(PT(1),PT(2),:,:,:,1));
%tmpINIT = shiftdim(squeeze(initPHI(PT(1),PT(2),:,:,:,1)),-1);
%waveT(:,:,2) = bsxfun(@plus,waveT(:,:,2),tmpINIT);

waveT(:,:,1) = 1;
%waveT(:,setdiff(1:size(waveT,2),4),:) = 0;

%waveT(:,4,2) = linspace(0,2*pi,size(waveT,1));


v_vec = ones(size(waveT));
sig_vec = pi*ones(size(waveT));

figure;
[wave] = buildWave(waveT,v_vec,sig_vec,theta,rho,R,NF,N);
figure;
imshow(oI(:,:,1),[]);
hold on
plot(PT(2),PT(1),'r*');
figure;
imshow(oI(:,:,1),[]);
%% look at dispersion measure and dynamcis
close all
for e = 1:size(retVec2,3)
    imshow(log(squeeze(retVec2(104,104,e,:,:,1))),[]);
    drawnow
waitforbuttonpress
end
figure;imshow(squeeze(retVec2(104,104,:,4,:,1)),[],'border','tight')
%%
SS = [];
for e = 1:13
    SS = cat(3,SS,GOUT2{1}{e});
end
%%
SSPP = [];
for e = 1:13
    SSPP = cat(3,SSPP,GOUT{e}{1});
end
%%
close all
h1 = figure;
for e = 1:size(SS,3)
    imshow([SS(:,130:150,e);GOUT2{1}(:,130:150,e)],[])

drawnow
pause(.3)
end
%% look at simple texture
for e = 1:size(oI,3)
    [d1 d2] = gradient(oI(:,:,e));
    d = (d1.^2 + d2.^2).^.5;
    TEX(e,:) = hist(d(:),linspace(-1,1,512));
end
%% gather histograms of gray scale images
for e = 1:size(oI,3)
    tmp = oI(:,:,e);
    [HI(e,:)] = hist(tmp(:),linspace(0,1,255));
end
%% cFlow over data
func = cFlow('gogo_bugEye');
func.setMCRversion('v930');
tsorFileList = issueBulkTicket(sorFileList);
serFea1 = {};
serFea2 = {};
serImg = {};
for e = 1:10
    [serFea1{e},serFeat2{e},serImg{e}] = func(tsorFileList{e},R,N,NF,mag,disp,func);
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,150,150);
%% decompose the feature 1
fprintf(['Start decomposing the whole.\n'])
fsel = 1;
numPC = 6;
sz = size(retVec1);
f1 = abs(squeeze(retVec1(:,:,:,:,1,:)));
clear retVec1
%f1 = diff(f1,1,3);
%f1 = retVec1;
sz1 = size(f1);
f1 = permute(f1,[1 2 3 5 4]);
%f1 = permute(f1,[1 2 4 5 3]);
sz11 = size(f1);
f1 = reshape(f1,[prod(sz11(1:4)) sz11(5)]);
[U,E,L] = PCA_FIT_FULLws(f1,numPC);
[C] = PCA_REPROJ(f1,E,U);
cI = reshape(C,[sz11(1:4) numPC]);
cI = permute(cI,[1 2 3 5 4]);
fprintf(['End decomposing the whole.\n'])
%% decompose the feature 2
fprintf(['Start decomposing the whole.\n'])
fsel = 1;
numPC = 6;
sz = size(retVec2);
%f1 = squeeze(retVec(:,:,:,:,:));
%f2 = abs(retVec2);
f2 = retVec2;
clear retVec2
%f2 = diff(retVec2,1,3);
sz1 = [size(f2) 1];
%f2 = reshape(f2,[sz1(1:3) prod(sz1(4:5)) sz1(6)]);
%f2 = reshape(f2,[sz1(1:3) prod(sz1(4)) sz1(5)]);
f2 = reshape(f2,[prod(sz1(1:2)) prod(sz1(3:4)) sz1(5)]);
%f2 = permute(f2,[1 2 4 5 3]);
%f2 = permute(f2,[1 2 3 5 4]);
f2 = permute(f2,[1 3 2]);
sz11 = size(f2);
f2 = reshape(f2,[prod(sz11(1:2)) sz11(3)]);
%f2 = zscore(f2);
[U2,E2,L2] = PCA_FIT_FULLws(f2,numPC);
[C2] = PCA_REPROJ(f2,E2,U2);
cI2 = reshape(C2,[sz11(1:4) numPC]);
cI2 = permute(cI2,[1 2 3 5 4]);
clear f2
%% decompose level 2 feature 2
numPC = 5;
f2 = single(cI2);
%f2 = diff(retVec2,1,3);
sz1 = [size(f2) 1];
f2 = permute(f2,[1 2 4 5 3]);
sz11 = size(f2);
f2 = reshape(f2,[prod(sz11(1:4)) sz11(5)]);
%f2 = zscore(f2);
[U3,E3,L3] = PCA_FIT_FULLws(f2,numPC);
[C3] = PCA_REPROJ(f2,E3,U3);
cIT2 = reshape(C3,[sz11(1:4) numPC]);
cIT2 = permute(cIT2,[1 2 3 5 4]);
%dK(isnan(dK(:)) | isinf(dK(:))) = 0;
%HUMM = PCA_REPROJ(dK(:)',E2,U2);
fprintf(['End decomposing the whole.\n']);
%% decompose level 2 feature 1
numPC = 5;
f2 = single(cI);
%f2 = diff(retVec2,1,3);
sz1 = [size(f2) 1];
f2 = permute(f2,[1 2 4 5 3]);
sz11 = size(f2);
f2 = reshape(f2,[prod(sz11(1:4)) sz11(5)]);
%f2 = zscore(f2);
[U4,E4,L4] = PCA_FIT_FULLws(f2,numPC);
[C4] = PCA_REPROJ(f2,E4,U4);
cIT1 = reshape(C4,[sz11(1:4) numPC]);
cIT1 = permute(cIT1,[1 2 3 5 4]);
%dK(isnan(dK(:)) | isinf(dK(:))) = 0;
%HUMM = PCA_REPROJ(dK(:)',E2,U2);
fprintf(['End decomposing the whole.\n']);
%% cluster the sampled feature space
f1 = cIT1;
sz = size(f1);
f1 = reshape(f1,[sz(1:2) prod(sz(3:4)) sz(5)]);
f1 = permute(f1,[1 2 4 3]);
sz = size(f1);
f1 = reshape(f1,[prod(sz(1:3)) sz(4)]);
f2 = cIT2;
sz = size(f2);
f2 = reshape(f2,[sz(1:2) prod(sz(3:4)) sz(5)]);
f2 = permute(f2,[1 2 4 3]);
sz = size(f2);
f2 = reshape(f2,[prod(sz(1:3)) sz(4)]);


fprintf(['Start clustering the part.\n'])
sam = 10;
clusterK = 7; % 7 for amplitude
usePC = [1 2 3 6 7 8 11 12];      % 3 for amplitude
options = statset('Display','iter','MaxIter',300);
GMModel = fitgmdist([f1(1:sam:end,usePC) f2(1:sam:end,usePC)],clusterK,'Options',options,'RegularizationValue',0.001,'Replicates',10);
%GMModel = fitgmdist([C2(1:sam:end,usePC)],clusterK,'Options',options,'RegularizationValue',0.00000001);
fprintf(['End clustering the part.\n'])
%% cluster the sampled feature space
fprintf(['Start clustering the part.\n'])
sam = 100;
clusterK = 5; % 7 for amplitude
usePC = [1:4];      % 3 for amplitude
options = statset('Display','iter','MaxIter',300);
GMModel = fitgmdist([C(1:sam:end,usePC) C2(1:sam:end,usePC)],clusterK,'Options',options,'RegularizationValue',0.001);
%GMModel = fitgmdist([C2(1:sam:end,usePC)],clusterK,'Options',options,'RegularizationValue',0.00000001);
fprintf(['End clustering the part.\n'])
%% cluster the whole data set via the gmm
fprintf(['Start clustering the whole via the part.\n'])
[IDX,~,PROB] = GMModel.cluster([f1(:,usePC) f2(:,usePC)]);
%[IDX,~,PROB] = GMModel.cluster([C(:,usePC) C2(:,usePC)]);
%[IDX,~,PROB] = GMModel.cluster([C2(:,usePC)]);

%IDXi = reshape(IDX,sz11(1:3));
IDXi = reshape(IDX,sz(1:3));
fprintf(['End clustering the whole via the part.\n']);
%% subcluster clicks
close all
cp = {};
S1 = reshape(cIT1,[size(cIT1,1) size(cIT1,2) size(cIT1,3)*size(cIT1,4) size(cIT1,5)]);
S2 = reshape(cIT2,[size(cIT2,1) size(cIT2,2) size(cIT2,3)*size(cIT2,4) size(cIT2,5)]);
St = cat(3,S1,S2);
SC = [];
toUse = [1 2 6 7 31 32 36 37];
toUse = [1 2 6 7 [1 2 6 7]+45];
%toUse = 1:60;
X = [];
Y = [];

NUMU = 40;
%%

h1 = figure;
h2 = figure;
toclick = false;
for e = 2:size(oI,3)
    
    
    
    
    
    SLICE = St(:,:,:,e);
    SLICE = reshape(SLICE,[size(SLICE,1)*size(SLICE,2) size(SLICE,3)]);
    if e ~= 1
        %{
        gm = gmdistribution(mean(SC,1),std(SC,1,1),1);
        [~,~,PROB] = gm.cluster(SLICE);
        PROB = reshape(PROB,[size(S1,1) size(S1,2)]);
        %}
        
        jSL = PCA_REPROJ(SLICE,jE,jU);
        [TYPE,PROB] = Mdl.predict(jSL);
        
        %{
        Values = SLICE(:,fIDX(1:30))*lambda;
        [TYPE,PROB] = Mdl.predict(Values);
        %}
        %{
        [TYPE,PROB] = tree.predict(SLICE(:,fIDX(1:3)));
        %}
        PROB = sum(PROB(:,1:cluster1),2);
        
        %{
        PROB = mvnpdf(SLICE(:,toUse),mean(SC(:,toUse)),cov(SC(:,toUse)));
        %{
        PROB = bsxfun(@minus,SLICE,mean(SC,1));
        PROB = bsxfun(@mtimes,PROB,std(SC,1,1).^-1);
        PROB = (sum(PROB.*PROB,2)).^.5;
        %}
        %}
        PROB = reshape(PROB,[size(S1,1) size(S1,2)]);
        PROB = log(.01*PROB);
        PROB(isinf(PROB)) = 0;
        PROB = imfilter(PROB,fspecial('gaussian',[31 31],5),'replicate');
        PROB = bindVec(PROB);
        
        
        if ~isempty(PROB)
            figure(h1);
            imshow(PROB,[]);
            if ~toclick
                figure(h2)
                out = flattenMaskOverlay(bindVec(oI(:,:,e)),reshape(TYPE,[size(oI,1) size(oI,2)])<= cluster1);
                imshow(out,[]);
                waitforbuttonpress
            end
            
        end
    end
  
    
    if toclick
        figure(h2)
        [cp{e}(:,2),cp{e}(:,1),~] = impixel(oI(:,:,e),[]);
        IDX = sub2ind([size(S1,1) size(S1,2)],cp{e}(:,1),cp{e}(:,2));
        msk = zeros(size(oI,1),size(oI,2));
        msk(IDX) = 1;
        msk = imdilate(msk,strel('disk',5,0));
        IDX = find(msk==1);
        Y = [Y;msk(:)];
        X = [X;SLICE];
        tS = [];
        for k = 1:size(St,3)
            tmp = St(:,:,k,e);
            tS = [tS tmp(IDX)];
        end

        SC = [SC;tS];


        if e >=1
            %[fIDX, fZ] = rankfeatures(X', Y,'Criterion','ttest');
            %Mdl = fitcnb(X(:,fIDX(1:NUMU)),Y);
            %tree = fitctree(X(:,fIDX(1:3)),Y,'PredictorSelection','curvature');
            %Mdl = fitcensemble(X(:,fIDX(1:3)),Y);
            
            
            %{
            [fIDX, fZ] = rankfeatures(X', Y,'Criterion','ttest');
            options = statset('Display','iter');
            
            [jU,jE,jL] = PCA_FIT_FULLws(X,50);
            jX = PCA_REPROJ(X,jE,jU);
            %}
            
            
            numF = 40;
            F1 = Y==1;
            
            AIC = [];
            BIC = [];
            GMModel_sub = {};
            for cluster1 = 1:10
                GMModel_sub{cluster1} = fitgmdist(jX(F1,:),cluster1,'Options',options,'Replicates',3);
                AIC(cluster1) = GMModel_sub{cluster1}.AIC;
                BIC(cluster1) = GMModel_sub{cluster1}.BIC;
            end
            [minAIC,cluster1] = min(AIC);
            cluster1 = 4;
        
            
           
            
            clear GMModel_sub;
            GMModel_sub = fitgmdist(jX(F1,:),cluster1,'Options',options,'Replicates',3);
            [subIDX1] = GMModel_sub.cluster(jX(Y==1,:));
            
            cluster2 = 5;
            F0 = find(Y==0);
            sam = 10;
            GMModel_sub = fitgmdist(jX(F0(1:sam:end),:),cluster2,'Options',options);
            [subIDX2] = GMModel_sub.cluster(jX(Y==0,:));
            
            YT = Y;
            YT(F1) = subIDX1;
            YT(F0) = subIDX2 + cluster1;
            Mdl = fitcnb(jX(:,:),YT);
            
            
            %{
            GMModel_sub = fitgmdist(X(Y==1,fIDX(1:10)),2,'Options',options);
            [subIDX] = GMModel_sub.cluster(X(Y==1,fIDX(1:10)));
            fidxY = find(Y==1);
            plot3(X(fidxY(subIDX==1),fIDX(1)),X(fidxY(subIDX==1),fIDX(2)),X(fidxY(subIDX==1),fIDX(3)),'r.')
            hold on
            plot3(X(fidxY(subIDX==2),fIDX(1)),X(fidxY(subIDX==2),fIDX(2)),X(fidxY(subIDX==2),fIDX(3)),'b.')
            
            GMModel_sub = fitgmdist(X(Y==0,fIDX(1:10)),2,'Options',options);
            [subIDX] = GMModel_sub.cluster(X(Y==0,fIDX(1:10)));
            fidxY = find(Y==0);
            plot3(X(fidxY(subIDX==1),fIDX(1)),X(fidxY(subIDX==1),fIDX(2)),X(fidxY(subIDX==1),fIDX(3)),'r.')
            hold on
            plot3(X(fidxY(subIDX==2),fIDX(1)),X(fidxY(subIDX==2),fIDX(2)),X(fidxY(subIDX==2),fIDX(3)),'b.')
            
            
            [lambda] = myLDA(X(:,fIDX(1:30)),Y);
            Values = X(:,fIDX(1:30))*lambda;
            Mdl = fitcnb(Values,Y);
            %}
        end
    end
end
%% view tensor oil
close all
h1 = figure;
h2 = figure;
for e = 1:size(IDXi,3)
    RGB = label2rgb(IDXin(:,:,e));
    figure(h1)
    imshow(RGB,[]);
    figure(h2)
    out = flattenMaskOverlay(bindVec(oI(:,:,e)),IDXin(:,:,e) == 8);
    imshow(out,[])
    waitforbuttonpress
end
%%
figure;
for e = 1:clusterK
    kidx = find(IDX == e);
    kidx = kidx(1:20:end);
    plot3(C2(kidx,1),C2(kidx,2),C2(kidx,3),'.');
    hold on
end
%% sub cluster after melt oil
fidx = find(IDX == 6);
clusterK = 3;
usePCS = 1:3;
usePC = [4 5 9 10];  
options = statset('Display','iter');
sam = 1;
%GMModel_sub = fitgmdist([C(fidx,usePCS) C2(fidx,usePCS)],clusterK,'Options',options);
GMModel_sub = fitgmdist([f1(1:sam:end,usePC) f2(1:sam:end,usePC)],clusterK,'Options',options);

[subIDX] = GMModel_sub.cluster([f1(fidx,usePC) f2(fidx,usePC)]);
nIDX = IDX;
nIDX(fidx) = subIDX + 7 -1;
IDXin = reshape(nIDX,sz(1:3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start over "time/scale" analysis
fprintf(['End clustering the whole via the part.\n'])
PROBi = reshape(PROB,[sz11(1:4) size(PROB,2)]);
PROBi = permute(PROBi,[1 2 5 4 3]);
pSZ = size(PROBi);
PROBi = reshape(PROBi,[prod(pSZ(1:4)) pSZ(5)]);
numPC_P = 3;
[pU,pE,pL] = PCA_FIT_FULLws(PROBi,numPC_P);
[pC] = PCA_REPROJ(PROBi,pE,pU);
fprintf(['End clustering the whole via the part.\n'])
%% try whole decomose - skip this for classic
pC = reshape(pC,[pSZ(1:4) size(pC,2)]);
pC = permute(pC,[1 2 3 5 4]);
hSZ = size(pC);
pC = reshape(pC,[hSZ(1:2) prod(hSZ(3:4)) hSZ(5)]);
%% try cluster via 2nd level decompose

tmpPC = permute(pC,[1 2 4 3]);
tmpSZ = size(tmpPC);
tmpPC = reshape(tmpPC,[prod(tmpSZ(1:3)) prod(tmpSZ(4))]);
fprintf(['Start clustering the part.\n'])
sam = 100;
clusterK = 3;
usePC = 3;
options = statset('Display','iter');
TMP_GMModel = fitgmdist(tmpPC(1:sam:end,1:usePC),clusterK,'Options',options);
fprintf(['End clustering the part.\n'])
kk = TMP_GMModel.cluster(tmpPC(:,1:usePC));
kk = reshape(kk,tmpSZ(1:3));
%% view 2nd level
close all
for e = 1:size(kk,3)
    RGB = label2rgb(kk(:,:,e));
    imshow(RGB,[]);
    drawnow
    waitforbuttonpress
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather the stats for the paths through the space
fprintf(['Start stats on the whole.\n'])
StateT = permute(IDXi,[1 2 4 3]);
stateSZ = size(StateT);
StateT = reshape(StateT,[prod(stateSZ(1:3)) stateSZ(4)]);
UQ = unique(IDX);


StateTT = im2col(StateT,[1 2],'sliding');
TT = [];
% gather trans probs
for u1 = 1:numel(UQ)
    for u2 = 1:numel(UQ)
        qS = [UQ(u1);UQ(u2)];
        O = all(bsxfun(@eq,StateTT,qS),1);
        TT(u1,u2) = sum(O) * size(StateTT,2)^-1;
        TT
    end
end

TT_r = [];
for r = 1:(size(StateT,2)-1)
    StateTT = im2col(StateT(:,r:r+1),[1 2],'sliding');
    % gather trans probs
    for u1 = 1:numel(UQ)
        for u2 = 1:numel(UQ)
            qS = [UQ(u1);UQ(u2)];
            O = all(bsxfun(@eq,StateTT,qS),1);
            TT_r(u1,u2,r) = sum(O) * size(StateTT,2)^-1;
        end
    end
end

% gather state probs
P = [];
UQ = unique(IDX);
for u = 1:numel(UQ)
    P(u) = sum(IDX==UQ(u)) * numel(IDX)^-1;
end
for u = 1:numel(UQ)
    rP(u,:) = sum(StateT == UQ(u),1) * size(StateT,1)^-1;
end
fprintf(['End stats on the part.\n'])
%% compute prob of patch
fprintf(['Start stats on the data.\n'])
STATEP_R = [];
for r = 1:size(StateT,2)
    ss = StateT(:,r);
    LOOKUP = rP(:,r);
    STATEP_R(:,r) = LOOKUP(ss)';
end

STATEP = [];
for r = 1:size(StateT,2)
    ss = StateT(:,r);
    STATEP(:,r) = P(ss)';
end

TRANSP = [];
for r = 1:(size(StateT,2)-1)
    ss = StateT(:,r:r+1);
    tIDX = sub2ind(size(TT),ss(:,1),ss(:,2));
    TRANSP(:,r) = TT(tIDX);
end


TRANS_rP = []
for r = 1:(size(StateT,2)-1)
    ss = StateT(:,r:r+1);
    LOOKUP = TT_r(:,:,r);
    tIDX = sub2ind(size(LOOKUP),ss(:,1),ss(:,2));
   
    TRANS_rP(:,r) = LOOKUP(tIDX);
end

totP = -(sum(log(STATEP),2) + sum(log(TRANSP),2));
totP_r = -(sum(log(STATEP_R),2) + sum(log(TRANS_rP),2));
probEfflux = -log(STATEP_R(:,1:end-1)) + -log(TRANS_rP);
fprintf(['End stats on the data.\n'])
%% plot state probs over radius
close all
figure;
CL = {'r' 'g' 'b' 'm' 'c' 'r--' 'g--'};
for e = 1:size(rP,1)
    plot(rP(e,:),CL{e})
    hold on
    LEG{e} = num2str(e);
end
legend(LEG)
%% decompose prob efflux
fprintf(['Start decompose on the efflux.\n'])
[oU,oE,oL] = PCA_FIT_FULLws(probEfflux,5);
[oilC] = PCA_REPROJ(probEfflux,oE,oU);
oilC = reshape(oilC,[sz11(1:2) NUM size(oilC,2)]);
fprintf(['End decompose on the efflux.\n'])
%% cluster oilC
oilRZ = size(oilC);
oC = reshape(oilC,[prod(oilRZ(1:3)) oilRZ(4)]);
fprintf(['Start clustering the part.\n'])
sam = 100;
clusterK = 4; % 7 for amplitude
options = statset('Display','iter','MaxIter',300);
GMModeOIL = fitgmdist([oC],clusterK,'Options',options,'RegularizationValue',0.00001);
kOIL = GMModeOIL.cluster(oC);
kOIL = reshape(kOIL,oilRZ(1:3));
%GMModel = fitgmdist([C2(1:sam:end,usePC)],clusterK,'Options',options,'RegularizationValue',0.00001);
fprintf(['End clustering the part.\n'])
%% reshape the prob for image
VW = reshape(totP,[sz11(1:2) sz11(4)]);
VW_r = reshape(totP_r,[sz11(1:2) sz11(4)]);
%% view prob of patch
s = 1
close all
figure;
imshow(bindVec(VW(:,:,s)),[]);

figure;
imshow(bindVec(VW_r(:,:,s)),[]);
figure;
imshow(cat(3,bindVec(oI(:,:,s)),zeros(size(oI,1),size(oI,2)),bindVec(VW(:,:,s))),[]);
figure;
imshow(oI(:,:,s),[]);
PP = bindVec(VW(:,:,s));
figure;
imshow(PP < .1)
%% view oil painting melt
close all
select = 1;
select = [4];
select = [2 3];
%select = [5 6 7 8];
%select = [2 9];
h1 = figure;
h2 = figure;
CL = {'r' 'b' 'g' 'm'};
for e = 1:size(IDXi,4)
    for r = 1:size(IDXi,3)
        tmp = IDXi(:,:,r,e);
        
       
        figure(h1)
        RGB = label2rgb(tmp);
        imshow(RGB,[],'Border','tight');
        title(num2str(e))
        drawnow
        
        figure(h2)
        out = oI(:,:,e);
        out = bindVec(out);
        for s = 1:numel(select)
            msk = tmp == select(s);
            out = flattenMaskOverlay(out,msk,.2,CL{s});
        end
        imshow(out,[],'Border','tight');
        
        if r ~= 1
            waitforbuttonpress
        end
            
        drawnow
    end
end
%% view abstract oil
close all
oilC= permute(oilC,[1 2 4 3]);
for e = 1:size(oilC,4)
    tmpO = oilC(:,:,:,e);
    for k = 1:size(tmpO,3)
        tmpO(:,:,k) = bindVec(tmpO(:,:,k));
    end
    imshow(tmpO,[]);
    figure;
    imshow(tmpO(:,:,1));
    figure;
    imshow(oI(:,:,e),[]);
    waitforbuttonpress
    close all
end
%% view raw image
close all
h1 = figure;
h2 = figure;
for e = 1:size(oI,3)
    figure(h1);
    imshow(oI(:,:,e),[]);
    figure(h2);
    RGB = label2rgb(kOIL(:,:,e));
    imshow(RGB,[]);
    drawnow
    
    waitforbuttonpress
end