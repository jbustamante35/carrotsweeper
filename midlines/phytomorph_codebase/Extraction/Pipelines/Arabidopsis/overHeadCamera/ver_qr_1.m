%% scan for NMS files - NEEDED FOR MAX_L
FilePath = '/mnt/tetra/nate/ABCD_Run_Two/';
FileList = {};
FileExt = {'jpg'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% sort scan for NMS files - NEEDED FOR MAX_L
FilePath = '/mnt/tetra/nate/ABCD_Run_Two/';
sFileList = {};
FileExt = {'jpg'};
sFileList = sdig(FilePath,sFileList,FileExt,1);
%% sample random images
close all
mag = .2;
FileList = FileList(randperm(numel(FileList)));
S = [];
cnt = 1;
e = 1;
while cnt < 100
    I = imread(FileList{e});
    I = imresize(I,mag);
    aI = mean(I(:));
    e = e + 1;
    if aI > 40
        sz = size(I);
        imshow(I,[]);
        drawnow
        I = reshape(I,[prod(sz(1:2)) sz(3)]);
        S = [S;I];
        cnt = cnt + 1;
        cnt
    end
end
S = double(S);
%% grab checkboards
close all
cS = [];
checkBoardList = {};
for e = 1:numel(sFileList)
    I = imread(sFileList{e}{1});
    checkBoardList{end+1} = sFileList{e}{1};
    I = imresize(I,mag);
    aI = mean(I(:));
    if aI > 40
        imshow(I,[]);
        title(sFileList{e}{1});
        drawnow
        sz = size(I);
        I = reshape(I,[prod(sz(1:2)) sz(3)]);
        cS = [cS;I];
        cnt = cnt + 1;
    end
end
cS = double(cS);
%% stack the random and checkerboard
mS = [S;cS];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% better fix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init phytoMorph typing graphs
mksqlite('close')
T = myDS();
imageFileNameType = T.putTypeNode('imageFileName');
fileNameType = T.putTypeNode('fileName');
referenceFrameType = T.putTypeNode('referenceFrame');
T.propLink(fileNameType,imageFileNameType,'typeof');
T.propLink(referenceFrameType,imageFileNameType,'hasa');
%%
nodeTable = table();
s = [];
for e = 1:50
    s(e).imageFileNameType = 'imageFileName';
    s(e).imageFileName = sFileList{50}{e};
    s(e).prop1 = 'referenceFrameType';
    s(e).link1 = 'hasa';
    s(e).data1 = eye(3);
end

%%
T.importImagesFromList(sFileList{50}(1:3),{'imageFileName','referenceFrameType',eye(3)},true);
%%
T.dijkstra(1,2,2)

%%
close all
for e = 1:numel(sFileList)
    toRead = 50;
    tf = 0;
    good = false;
    while ~good
        tFile{e} = sFileList{e}{toRead};
        I = imread(tFile{e});
        tf = isNight(I,50);
        if ~tf
            imshow(I);
            good = true;
        else
            toRead = toRead + 10;
            good = false;
        end
        toRead
    end
end
%% collect data for zoom nets
close all
toUse = 10;
pointList = zoomGather_nonDB(tFile(1:toUse),[200 200],true);
%% collect data for sparkel nets
[sparklePointList,sparkleUse] = zoomGather_nonDB(tFile(1:50),[200 200],false);
%% sample zoom sequence
boxSequence = [[250 250];[150 150]];
[X,Y] = sampleZoomSequence(tFile(1:toUse),pointList,boxSequence,[.25 .75],[100 50],25,[-pi/32,pi/32,100],false);
%% sample sparkle
[sparkleX,sparkleY] = sampleSparkleSequence(tFile(1:50),.05,sparklePointList,true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train zoom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nHood = [];
nHood(1,:) = [4 4];
nHood(2,:) = [4 4];
for z = 1:numel(X)
    tX = X{z};
    tY = Y{z};
    close all
    %%%%%%%%%%%%%%%%%%%%%
    layers = [imageInputLayer([size(tX,1) size(tX,2) size(tX,3)],'Normalization','None');
              convolution2dLayer(nHood(z,:),11);
              reluLayer();
              maxPooling2dLayer(2,'Stride',2);
              fullyConnectedLayer(size(tY,1));
              regressionLayer();];
    options = trainingOptions('sgdm',...
        'InitialLearnRate',.005,...
        'MaxEpochs',2000,...
        'Plots','training-progress',...
        'ExecutionEnvironment','parallel');
    [tY,uY{z},sY{z}] = zscore(tY');
    %tY = tY';
    upperRightCorner{z} = trainNetwork(tX,tY,layers,options);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% meta train zoom - local
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for z = 1:numel(X)
    
    tX = X{z};
    tY = Y{z};
    close all
   
    %uY{z} = mean(tY',1);
    %tY = bsxfun(@minus,tY',uY{z});
    [tY,uY{z},sY{z}] = zscore(tY');

    N = 750;
    IDX = randperm(size(tX,4));
    trainX = tX(:,:,:,IDX(1:N));
    trainY = tY(IDX(1:N),:);
    testX = tX(:,:,:,IDX((N+1):end));
    testY = tY(IDX((N+1):end),:);


    nHood = optimizableVariable('nHood',[3,15],'Type','integer');
    nKernels = optimizableVariable('nKernels',[3,15],'Type','integer');
    initLearnRate = optimizableVariable('initLearnRate',[.0001,.1]);
    L2Regularization = optimizableVariable('L2Regularization',[.0001,.0008]);
    Momentum = optimizableVariable('Momentum',[0,1]);
    exeEnvironment = 'parallel';
    
    para = [nHood,nKernels,initLearnRate,L2Regularization,Momentum];
    
    [BO2{z},net{z}] = gpuTrainRegression(trainX,trainY,testX,testY,exeEnvironment,para);
   
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% meta train zoom - remote
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nTrain = [300 5000];
for z = 1:numel(X)
    
    func = cFlow('gpuTrainRegression');
    func.setMCRversion('v930');
    func.setMemory('8000');
    func.setGPU(1);
    
    
    tX = X{z};
    tY = Y{z};
    close all;
   
    %uY{z} = mean(tY',1);
    %tY = bsxfun(@minus,tY',uY{z});
    [tY,uY{z},sY{z}] = zscore(tY');

    N = round(.75*size(tX,4));
    
    
    IDX = randperm(size(tX,4));
    trainX = tX(:,:,:,IDX(1:N));
    trainY = tY(IDX(1:N),:);
    testX = tX(:,:,:,IDX((N+1):end));
    testY = tY(IDX((N+1):end),:);


    nHood = optimizableVariable('nHood',[7,15],'Type','integer');
    nKernels = optimizableVariable('nKernels',[3,15],'Type','integer');
    initLearnRate = optimizableVariable('initLearnRate',[.0001,.1]);
    L2Regularization = optimizableVariable('L2Regularization',[.0001,.0008]);
    Momentum = optimizableVariable('Momentum',[0,1]);
    exeEnvironment = 'gpu';
    
    para = [nHood,nKernels,initLearnRate,L2Regularization,Momentum];
    
    
    [BO{z},net{z}] = func(nTrain,trainX,trainY,testX,testY,exeEnvironment,para);
    
    
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    func.submitDag(auth,50,50);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loader
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for z = 1:numel(X)
    upperRightCorner{z} = cFlowLoader(net{z});
end
%% test loaded on data
z = 1;
for z = 1
    for e = 1:100
        tX = X{z};
        pY = upperRightCorner{z}.predict(tX(:,:,:,e));
        pY = uY{z} + sY{z}.*pY;
        tY = Y{z};
        imshow(tX(:,:,:,e),[]);
        hold on
        drawnow

        plot(tY(1,e),tY(2,e),'g.');
        plot(pY(1),pY(2),'go');
        drawnow
        waitforbuttonpress
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train sparkle - single local
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tY = sparkleY;
tX = sparkleX;
layers = [imageInputLayer([size(tX,1) size(tX,2) size(tX,3)],'Normalization','None');
          convolution2dLayer([11 11],10);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(size(tY,2));
          regressionLayer();];
options = trainingOptions('sgdm',...
    'InitialLearnRate',.001,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');
%uY{z} = mean(tY',1);
%tY = bsxfun(@minus,tY',uY{z});
[tY,suY,ssY] = zscore(tY);
%tY = tY';
usparkle = trainNetwork(tX,tY,layers,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% meta train sparkle - remote
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nTrain = [300 5000];

func = cFlow('gpuTrainRegression');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);

tY = sparkleY;
tX = sparkleX;
close all

[tY,suY,ssY] = zscore(tY);

N = 40;

IDX = randperm(size(tX,4));
trainX = tX(:,:,:,IDX(1:N));
trainY = tY(IDX(1:N),:);
testX = tX(:,:,:,IDX((N+1):end));
testY = tY(IDX((N+1):end),:);

nHood = optimizableVariable('nHood',[7,15],'Type','integer');
nKernels = optimizableVariable('nKernels',[3,15],'Type','integer');
initLearnRate = optimizableVariable('initLearnRate',[.0001,.1]);
L2Regularization = optimizableVariable('L2Regularization',[.0001,.0008]);
Momentum = optimizableVariable('Momentum',[0,1]);
exeEnvironment = 'gpu';

para = [nHood,nKernels,initLearnRate,L2Regularization,Momentum];

[BO_sparkel,net_sparkle] = func(nTrain,trainX,trainY,testX,testY,exeEnvironment,para);

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%%
usparkle = cFlowLoader(net_sparkle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test sparkle zoom
%close all
I = double(imread(tFile{end-1}))/255;
sI = imresize(I,.05);

pY = [];
pY = fitbitp(sI,4,Smap,Sbeta,SWINDOW);
%{
pY = usparkle.predict(sI);
pY = bsxfun(@times,pY,ssY);
pY = bsxfun(@plus,pY,suY);
%}
pY = reshape(pY,[9 2]);

zoomSequence = [.25 .75];
ITER = [5 5];
for p = 1:size(pY,1)
    cnt = 1;
    for z = 1:2
        for iter = 1:ITER(z)
            
            
            BOX = point2Box((double(round(pY(p,:,cnt)))),boxSequence(z,:));
            [subI] = mimcrop(I,BOX,zoomSequence(z),[]);
            
            
            POS = fitbitp(subI,4,map{z},beta{z},WINDOW{z});
            deltaZ = POS - [size(subI,2) size(subI,1)]/2;
            
            %deltaZ = uY{z} + sY{z}.*upperRightCorner{z}.predict(subI) - [size(subI,2) size(subI,1)]/2;
            pY(p,:,cnt+1) = pY(p,:,cnt) + zoomSequence(z)^-1*((deltaZ));
            cnt = cnt + 1;
        end
    end
end
figure;
imshow(I,[]);
hold on

for p = 1:size(pY,1)
    plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'g')
end
plot(pY(:,1,1),pY(:,2,1),'g.');
plot(pY(:,1,end),pY(:,2,end),'m*');
%% fitbit LOCK 
clear map beta WINDOW
for z = 1:numel(X)
    [map{z},beta{z},WINDOW{z}] = fitbit(X{z},Y{z}',4,64);
end
% fitbit SPARKLE 
clear Smap Sbeta SWINDOW
[Smap,Sbeta,SWINDOW] = fitbit(sparkleX,sparkleY,4,64);

%% dither test

[dX,map] = reduce(X{1},4);
RGB = ind2rgb(dX(:,:,4),map);
[bX] = bits(dX,4);
%% 
[map,beta,WINDOW] = fitbit(X{1},Y{1}',4,64);
%% try fitbitp
Y = fitbitp(X{1}(:,:,:,1:3),4,map,beta,WINDOW);
%%
I = imread(sFileList{50}{3});
%% stack wrap
HSV = rgb2hsv(mS/255);
data = HSV(:,1);
loc = 1;
mx = 1;
mdata2 = data+mx;
mdata1 = data-mx;
dd = [(mdata1);(data);(mdata2)];
%dd = dd(dd > -.2 & dd < .17);
close all
hist(dd,linspace(-1,2,20000));
figure;
hist(mS(:,2),linspace(0,255,256));
%% unmix the last channel
GMModel = fitgmdist(mS(:,2),4);
%%
close all
subL = GMModel.cluster(mS(:,2));
idx = subL==4 & HSV(:,1) > .15 & HSV(:,1) < .3;% & HSV(:,2) > .3 & HSV(:,2) < .8 ;
hist(HSV(idx,2),linspace(0,1,255))
%%
GMModel = fitgmdist(HSV(subL==4,4));
%% try to label value channel for an image
close all
I = imread(sFileList{e}{1});
sz = size(I);
figure;
imshow(I,[]);
HSVi = rgb2hsv(double(I)/255);
figure
imshow((HSVi(:,:,1) < .05 | HSVi(:,:,1) > .97) & HSVi(:,:,3) > .4 ,[]);
figure;
imshow((HSVi(:,:,1) > .15 & HSVi(:,:,1) < .3) & HSVi(:,:,2) > .189,[]);


figure;
imshow((I(:,:,2) > 40 & I(:,:,2) < 255) & (HSVi(:,:,1) > .15 & HSVi(:,:,1) < .3),[]);
%%
figure;
imshow(HSVi(:,:,2) < .3 & HSVi(:,:,3) > .8,[]);
figure;
imshow((HSVi(:,:,1) > .3 & HSVi(:,:,1) < .7) & HSVi(:,:,2) > .38 & HSVi(:,:,3) < .77,[]);
figure;
imshow(HSVi,[]);
figure;
imshow(HSVi(:,:,3),[])
figure;
imshow(HSVi(:,:,2),[])
figure;
imshow(HSVi(:,:,1),[])
figure;
imhist(HSV(:,2));
figure;
imhist(HSV(:,3));
figure;
imhist(HSV(:,1));
%%
I = reshape(I,[prod(sz(1:2)) sz(3)]);
HSVi = rgb2hsv(double(I)/255);
cidx = GMModel.cluster(HSVi(:,3));
cidx = reshape(cidx,sz(1:2));
close all
imshow(cidx,[])
%% cluster the rgb values
clusterK = 8;
sub = 100;
hsv = rgb2hsv(mS(1:sub:end,:)/255);
[dis] = myWrapperDistance(hsv(:,1),.9,1);
func = @(X)myWrapperG(hsv(:,1),X);
ops = optimset('Display','iter');
ops = optimoptions('fmincon','Display','iter');
x = fmincon(func,[.1 .25 .55 .3 .3 .3],[],[],[],[],[0 0 0 0 0 0],[1 1 1 100 100 100],[],ops);
x = fminsearch(func,[.1 .25 .55 .3 .3 .3],ops);

MY = @(X,Y)myPF1(hsv(:,1),X,Y);
ops = optimoptions('particleswarm','Display','iter','UseParallel',true,'PlotFcn',MY);
x = particleswarm(func,6,[0 0 0 0 0 0],[1 1 1 .3 .3 .3],ops);


[phat,pci] = mle(x,'pdf',@(x,v,d)ncx2pdf(x,v,d),'start',[1,1])
options = statset('Display','iter','MaxIter',300);
GMModel = fitgmdist([mS(1:sub:end,:)],clusterK,'Options',options,'RegularizationValue',0.00001,'Replicates',3);
%% stack random list and checkboard list
rList = {};
for e = 1:numel(checkBoardList)
    rList{end+1} = checkBoardList{e};
end
for e = 1:500
    rList{end+1} = FileList{e};
end
%%
close all force
msg{1} = 'Please click on red';
msg{2} = 'Please click on plant';
msg{3} = 'Please click on black background';
msg{4} = 'Please click on white QR code';
msg{5} = 'Please click on blue checkboard border';
I = imread(sFileList{70}{1});
figure;
imshow(I,[])

figure;

L = labelImage(double(I),GMModel);
imshow(L,[]);
%V = {[5] [7] [8 1 4] [2] [3]};
V = {[1] [2] [5] [3] [4]};
%%
V = [];
for e = 1:numel(msg)
    h = msgbox(msg{e});
    [~,~,V(e,:)] = impixel(uint8(L),[]);
end
V = mean(V,2);
CL = {'r' 'g' 'k' 'w' 'b'};
close all force
%%
disp = false;
close all
oPath = '/mnt/tetra/nate/ABCD_Run_Two/return3/';
oPath = '/mnt/tetra/nate/overHeadReturn/';
mkdir(oPath);
CL = {'r' 'g' 'k' 'w' 'b'};
for s = 1:numel(sFileList)
    try
        
        [C] = getPlantCells(sFileList{s},[],V,CL,oPath,disp);
        
    catch ME
        ME
    end
end
%% whole run
FilePath = '//mnt/tetra/nate/overHead/';
sFileList = {};
FileExt = {'jpg'};
sFileList = sdig(FilePath,sFileList,FileExt,1);
%% look at first
% cook
%%
close all
imshow(nI,[]);
%%
close all

RGB = label2rgb(L);
imshow(RGB,[]);
%% show movie
s  =1;
for e = 1:numel(sFileList{s})
    I = imread(sFileList{s}{e});
    imshow(I,[]);
drawnow
end
%% apply to random
close all
for e = 1:numel(rList)
    oI = imread(rList{e});
    if mean(double(oI(:))) > 40
        sz = size(oI);
        I = reshape(oI,[prod(sz(1:2)) sz(3)]);
        I = double(I);
        label = GMModel.cluster(I);
        label = reshape(label,sz(1:2));
        
        
        
        
        redTape = label == 1;
        plant = label == 5;
        blackBackground = label == 3 | label == 4;
        whiteBackground = label == 2;
        checkBoard = label == 6;
        
        hcheckBoard = imfill(checkBoard,'holes');
        ccheckBoard = hcheckBoard == 1  & checkBoard == 0;
        ccheckBoard = bwlarge(ccheckBoard);
        checkA = sum(ccheckBoard(:));
        
        if checkA > 500000
            checkerFlag = true;
            out = flattenMaskOverlay(double(oI)/255,ccheckBoard,.3,'b');
            out = flattenMaskOverlay(out,checkBoard,.8,'b');
            
            [imagePoints,boardSize] = detectCheckerboardPoints(double(oI)/255);
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % order the checkerboard points
            kidx = kmeans(imagePoints(:,2),5);
            L = [];
            for g = 1:5
                L = [L ,kidx==g];
                v(g) = mean(imagePoints(kidx==g,2));
            end
            [v,sidx] = sort(v);
            L = L(:,sidx);
            resEstimate = mean(diff(v,1,2));
            % my sort
            iS = [];
            for g = 1:5
                fidx = L(:,g)==1;
                d = imagePoints(fidx,:);
                [~,sidx] = sort(d(:,1));
                d = d(sidx,:);
                iS = [iS;d];
            end
            imagePoints = iS;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            worldPoints = generateCheckerboardPoints([6 6],1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % order the checkerboard points
            kidx = kmeans(worldPoints(:,2),5);
            L = [];
            for g = 1:5
                L = [L ,kidx==g];
                v(g) = mean(worldPoints(kidx==g,2));
            end
            [~,sidx] = sort(v);
            L = L(:,sidx);

            % my sort
            iS = [];
            for g = 1:5
                fidx = L(:,g)==1;
                d = worldPoints(fidx,:);
                [~,sidx] = sort(d(:,1));
                d = d(sidx,:);
                iS = [iS;d];
            end
            worldPoints = iS;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            worldPoints = bsxfun(@minus,worldPoints,mean(worldPoints,1));
            %worldPoints = bsxfun(@plus,worldPoints,[size(oI,2) size(oI,1)]);
            
            
            params = estimateFisheyeParameters(cat(3,imagePoints,imagePoints),worldPoints,[size(oI,1) size(oI,2)]);
            J = undistortFisheyeImage(oI,params.Intrinsics);
            
            
            tform = fitgeotrans(worldPoints,imagePoints,'projective');
            
            
            W = round(.5*(size(oI,1).*resEstimate^-1));
            NIP = round(2*W*resEstimate);
            [w1 w2] = ndgrid(linspace(-W,W,NIP),linspace(-W,W,NIP));
            XI = tform.transformPointsForward([w2(:),w1(:)]);
            nI = [];
            for k = 1:size(oI,3)
                nI(:,k) = ba_interp2(double(oI(:,:,k)),XI(:,1),XI(:,2));
            end
            nI = reshape(nI,[NIP NIP 3]);
            
            
            figure;
            imshow(nI/255,[]);
            figure
            imshow(oI);
            waitforbuttonpress
         
            
            
            
        else
            checkerFlag = false;
            
            redTape = imclose(redTape,strel('line',400,90));
            redTape = imclose(redTape,strel('line',400,0));


            out = flattenMaskOverlay(double(oI)/255,redTape,.3,'r');
            out = flattenMaskOverlay(out,plant,.3,'g');
            out = flattenMaskOverlay(out,blackBackground,.3,'k');
            out = flattenMaskOverlay(out,whiteBackground,.3,'w');
            out = flattenMaskOverlay(out,checkBoard,.3,'b');



            innerCell = imfill(redTape,'holes');
            innerCell = innerCell == 1 & redTape == 0;
            innerCell = bwlarge(innerCell,9);
            R = regionprops(innerCell,'BoundingBox');


            RGB = label2rgb(label);
            
            
        end
        
        imshow(out,[]);
        hold on
        
        % display section
        if ~checkerFlag
            for b = 1:numel(R)
                rectangle('Position',R(b).BoundingBox,'EdgeColor','b','LineWidth',3);
            end
        else
            plot(imagePoints(:,1),imagePoints(:,2),'b*')
            plot(imagePoints(:,1),imagePoints(:,2),'ko')
        end
        drawnow
        hold off
        %waitforbuttonpress
    end
end
