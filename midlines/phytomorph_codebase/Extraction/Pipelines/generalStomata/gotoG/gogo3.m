%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the data for each species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issue ticket(s) - issue write ticket(s) over the directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TmaizeFileList = issueBulkTicket(maizeFileList);
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST_CORN2/';
[rPat, iticket] = issueTicket(rPath(1:end-1),10*numel(TmaizeFileList),'write');
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
TsetFileList = issueBulkTicket(setFileList);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the data for each species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather the normalization data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%massDownload('/iplant/home/phytomorphuser/workIT/', '.mat','/mnt/tetra/nate/workIT/');
massDownload('/iplant/home/phytomorphuser/workITOUT_FAST_CORN2/', '.mat','/mnt/tetra/nate/workIT_CORN_TRAIN/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scan local for mat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cMatFilePath = '/mnt/tetra/nate/workIT/';
cMatFilePath = '/mnt/tetra/nate/workIT_CORN_TRAIN/';
cMatFilePath = '/mnt/tetra/nate/workIT_SET_TRAIN/';
fPatchFileList = {};
FileExt = {'mat'};
fPatchFileList = gdig(cMatFilePath,fPatchFileList,FileExt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build out basis vector - from CONDOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[basisU,basisE] = buildGradeBasis(fPatchFileList,100);
%[basisU_CORN,basisE_CORN] = buildGradeBasis(fPatchFileList,100);
[basisU_SET,basisE_SET] = buildGradeBasis(fPatchFileList,100);
%[data] = myRemoteLoader(gradeFileList{1},'T');
%data = prepareData(fPatchFileList{1},basisU,basisE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build out human click map(s) - from CONDOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a temp local location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpPath = ['/mnt/tetra/nate/tmp/' tempname filesep];
mkdir(tmpPath);
humanClickPath = '/mnt/tetra/nate/sorStomataClick/';
humanClickPath = '/mnt/tetra/nate/maizeStomataClick/';
humanClickPath = '/mnt/tetra/nate/maizeStomataClick/';
mkdir(humanClickPath);
toClick = 100;
for e = 1:toClick
    try
        tmpM = zeros(size(rawI));
        [pth,nm,ext] = fileparts(fPatchFileList{e});
        outFile = [humanClickPath nm '.mat'];
        fileName = fPatchFileList{e};
        if ~exist(outFile)
            fileList{1} = fileName;
            %fileList = xfer_get({fileName},tmpPath,0,0);
            % change from rI to rawI at UIUC - not sure why i renamed the
            % varaiable-  but i did
            load(fPatchFileList{e},'T','rawI','customData');

            tmpM(customData{2}) = 1;
            R = regionprops(logical(tmpM),'boundingBox');
            subI = imcrop(rawI,R(1).BoundingBox);
            uix = [];
            [uix(:,2),uix(:,1),~] = impixel(subI,[]);

            
            rI = subI;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build out X and Y for training - from CONDOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
humanClickPath = '/mnt/tetra/nate/sorStomataClick/';
humanClickPath = '/mnt/tetra/nate/maizeStomataClick/';
humanClickPath = '/mnt/tetra/nate/seteriaStomataClick/';
trainFileList = {};
FileExt = {'mat'};
trainFileList = gdig(humanClickPath,trainFileList,FileExt,1);
dilateValue = 9;
%[labeledTrainingPackage] = generateMaskPackage(trainFileList,100,dilateValue,basisU,basisE);
[labeledTrainingPackage_CORN] = generateMaskPackage(trainFileList,100,dilateValue,basisU_CORN,basisE_CORN);
%% train labeled package - CONDOR workflow
%[AI_layer] = trainAIlayer(labeledTrainingPackage,[1 1]);
[AI_layer_CORN] = trainAIlayer(labeledTrainingPackage_CORN,[1 1]);
%% applyAI layer - CONDOR
[probMap] = applyAIlayer(fPatchFileList{1},AI_layer,basisU,basisE,[108 108]);
%% optimize AI output
opti_para = {};
%% optimize output layer with swarm and nelder mead
%[opti_para] = optimizeAIoutput(AI_layer,labeledTrainingPackage);
[opti_para_CORN] = optimizeAIoutput(AI_layer_CORN,labeledTrainingPackage_CORN,opti_para_CORN);
%% test opti
[probMap] = applyAIlayer(fPatchFileList{4},AI_layer,basisU,basisE,[108 108]);
[grade,ret1,ret2,BO] = opti(probMap,'',opti_para,labeledTrainingPackage.oI(1:108,1:108,4));
%%
applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,{'hello'},AI_layer,basisU,basisE,opti_para,'','');
%%
func = @(X)applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,{'hello'},AI_layer,basisU,basisE,opti_para,'./output/','');
applyFuncWrapper = partialFunction(func,'metaLayers_sorghum');
%applyFuncWrapper.publish();
%%
func_CORN = @(X)applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU_CORN,basisE_CORN,AI_layer_CORN},domainData,pointSet,{'hello'},AI_layer_CORN,basisU_CORN,basisE_CORN,opti_para_CORN,'./output/','');
applyFuncWrapper = partialFunction(func_CORN,'metaLayers_corn');
applyFuncWrapper.publish();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate training data on CONDOR for sorghum data
% need to run through feature extraction first before clicking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate operational grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I = imread(sorFileList{1});
%I = imread(maizeFileList{1});
I = imread(setFileList{1});
imageSZ = size(I);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issue ticket(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TsorFileList = issueBulkTicket(sorFileList);
%TmaizeFileList = issueBulkTicket(maizeFileList);
TsetFileList = issueBulkTicket(setFileList);
rPath = '/iplant/home/nmiller/workIT/';
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST/';
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST_CORN2/';
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST_SET/';
[rPath,iticket] = issueTicket(rPath(1:end-1),10*numel(setFileList),'write');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% issue ticket(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = './outputFeatures/';
fibratedExtrationLayer = cFlow('extractionLayer');
fibratedExtrationLayer.setMCRversion('v930');
for e = 1:300
    for sD = 1%:numel(subPointSet)
        %[pth,nm,ext] = fileparts(TsorFileList{e});
        %[pth,nm,ext] = fileparts(TmaizeFileList{e});
        [pth,nm,ext] = fileparts(TsetFileList{e});
        customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(sD) '}'];
        customData{2} = globalIDX{sD};
        % for local run
        %fibratedExtrationLayer(TsorFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
        %fibratedExtrationLayer(TmaizeFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
        fibratedExtrationLayer(TsetFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
    end
    e
end
fibratedExtrationLayer.submitDag(auth,150,150);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run all layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = './outputFeatures/';
func = cFlow('cRunner');
func.setMCRversion('v930');
for e = 1:100%numel(TmaizeFileList)
    for sD = 1%:numel(subPointSet)
        % set file name
        [pth,nm,ext] = fileparts(maizeFileList{e});
        customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(sD) '}'];
        %customData{2} = globalIDX{sD};
        % for local run
        func(maizeFileList{e},oPath,rPath,customData);
       
        sD
    end
    e
end
func.submitDag(auth,400,400);
%%
%%%%%%%%%%%
% to publish generalized loader
applyFunc = @(X,oPath,rPath,OTHER)applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,OTHER,AI_layer,basisU,basisE,opti_para,oPath,rPath);
applyFuncWrapper = partialFunction(applyFunc,'sorghumStomataApply');
applyFuncWrapper.publish();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% single application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 1;
rtPath = '';
otPath = '';
[fdata] = applyAllLayers(sorFileList{e},'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,customData,AI_layer,basisU,basisE,opti_para,otPath,rtPath);
%%
cRunner(TsorFileList{e},oPath,rPath,customData);