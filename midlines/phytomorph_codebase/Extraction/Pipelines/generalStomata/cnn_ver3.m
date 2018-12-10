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
%% NEW STEP:0
I = imread(iFileList{1});
cropSZ = 40/2;
%subP = {};
close all
R = 40;
for e = 3%2%:numel(maizeStomataCenter_row)
   
    tmpI = imread(iFileList{e});
    tmp = zeros(size(tmpI));
    imshow(tmpI,[]);
    hold on
    plot(maizeStomataCenter_row{e},maizeStomataCenter_column{e},'.')
    [q1 q2 V] = impixel();
    new{e} = [[maizeStomataCenter_row{e};q1],[maizeStomataCenter_column{e};q2]];
    ridx = find(any(new{e}(:,1) < R | new{e}(:,2) < R | new{e}(:,1) > (size(tmpI,1) - R) | new{e}(:,2) > (size(tmpI,2) - R),2));
    bad{e} = new{e}(ridx,:);
    new{e}(ridx,:) = [];
    plot(new{e}(:,1),new{e}(:,2),'ro')
    plot(bad{e}(:,1),bad{e}(:,2),'g*')
    waitforbuttonpress
    close all
  
    for p = 1:size(new{e},1)
         BOX = [new{e}(p,1)-cropSZ new{e}(p,2)-cropSZ cropSZ*2 cropSZ*2];
         sub = imcrop(tmpI,BOX);
         sub = interp2(sub,4);
         [subP{e}{p}(:,1) subP{e}{p}(:,2) V] = impixel(sub);
         subP{e}{p} = subP{e}{p} * .25;
    end
end
%% NEW STEP:1
e=1;
p=1;
e = 1;

SAMP = [];
cnt = 1;
storePARA = {};
rot = linspace(-pi/2,pi/2,100);
s1 = linspace(.8,1.2,100);
s2 = linspace(.8,1.2,100);
selI = [];
cntR = 1;
cntS1 = 1;
cntS2 = 1;
rotI = [];
yROT = [];
for e = 1:3
    I = imread(iFileList{e});
    iP = sub2ind(size(I),new{e}(:,2),new{e}(:,1));
    for p = 1:numel(subP{e})
        %close all
       
        [PARA] = generateTransformationPara(subP{e}{p},cropSZ);
        %PARA(4:5) = new{e}(p,:);
        storePARA{e}(p,:) = PARA;
        [T] = generateDomainTransformation(PARA);

        [t1 t2] = ndgrid(linspace(-2,2,200),linspace(-2,2,200));
        tD = [t2(:) t1(:)];

        [t1 t2] = ndgrid(linspace(0,2,100),linspace(-pi,pi,400));
        tD = [t1(:).*cos(t2(:)) t1(:).*sin(t2(:))];

        rsz = size(I);
        [tD] = transformDomain(tD,T);
        [sI] = myInterp2Sampler(I,iP(p)',tD,size(t1));


        % sample the square grid
        [x1 x2] = ndgrid(linspace(-25,25,51),linspace(-25,25,51));
        xD = [x2(:) x1(:)];

        % get many angle for guessing
        for r = 1:numel(rot)
            xD = [x2(:) x1(:)];
            tmpPARA = PARA;
            tmpPARA(1) = tmpPARA(1) + rot(r);
            tmpPARA(2:3) = 1;
            [T] = generateDomainTransformation(tmpPARA);
            [xD] = transformDomain(xD,T);
            [xI] = myInterp2Sampler(I,iP(p)',xD,size(x1));
            szr = size(xI);
            xI = zscore(xI(:));
            xI = reshape(xI,szr);
            
            if r == 1
                %imshow(xI,[]);
                %drawnow
            end
            rotI(:,:,cntR) = xI;
            %waitforbuttonpress
            cntR = cntR + 1;
            %tmpPARA

        end
        
        
        % get many stretch-1
        for r = 1:numel(s1)
            xD = [x2(:) x1(:)];
            tmpPARA = PARA;
            yS1(cntS1) = tmpPARA(2)*s1(r)^-1;
            tmpPARA(2:3) = 1;
            tmpPARA(2) = s1(r);
            [T] = generateDomainTransformation(tmpPARA);
            [xD] = transformDomain(xD,T);
            [xI] = myInterp2Sampler(I,iP(p)',xD,size(x1));
            %imshow(xI,[]);
            %drawnow
            SAM1(:,:,cntS1) = xI;
            
            cntS1 = cntS1 + 1;
            %tmpPARA
            
        end
        
         % get many stretch-2
        for r = 1:numel(s2)
            xD = [x2(:) x1(:)];
            tmpPARA = PARA;
            yS2(cntS2) = tmpPARA(3)*s2(r)^-1;
            tmpPARA(2:3) = 1;
            tmpPARA(3) = s2(r);
            [T] = generateDomainTransformation(tmpPARA);
            [xD] = transformDomain(xD,T);
            [xI] = myInterp2Sampler(I,iP(p)',xD,size(x1));
            %imshow(xI,[]);
            %drawnow
            SAM2(:,:,cntS2) = xI;
            
            cntS2 = cntS2 + 1;
            %tmpPARA
            
        end
        
        
        
        
        yROT = [yROT;rot'];
        brot = zeros(size(rot));
        brot(1) = 1;
        selI = [selI;brot'];


        TH = linspace(-pi,pi,100);
        p1 = 50*cos(TH);
        p2 = 50*sin(TH);

        %{
        figure;
        imshow(I,[]);
        hold on
        plot(new{e}(p,1),new{e}(p,2),'r*');

        drawnow
        figure;
        imshow(sI,[]);
        hold on
        plot(p1+size(sI,1)/2,p2+size(sI,2)/2,'r');
        title(num2str(PARA(1)*180/pi))
        pause(.001);
        drawnow
        %}

        SAMP(:,:,cnt) = sI;
        cnt = cnt + 1;
        %waitforbuttonpress


    end
end
%%
yS2(end) = [];
yS1(end) = [];
%%
close all
for e = 1:100
    imshow(rotI(:,:,e),[]);
    drawnow
    pause(.2)
end
%%
close all
for e = 1:100
    imshow(SAM1(:,:,e),[]);
    drawnow
    pause(.1)
end
%%
close all
for e = 1:100
    imshow(SAM2(:,:,e),[]);
    drawnow
    pause(.1)
end
%%
rr = [];
for e = 1:size(rotI,3)
    rr(:,:,e) = imresize(rotI(:,:,e),.25);
    e
end
%%

%% regress angle
rr = rotI;
layers = [ ...
    imageInputLayer([size(rr,1) size(rr,2) 1])
    convolution2dLayer([11 11],7)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','InitialLearnRate',0.01, 'MaxEpochs',150,'Plots','training-progress','ExecutionEnvironment','parallel');
xROT = reshape(rr,[size(rr,1) size(rr,2) 1 size(rr,3)]);
[newY,zu0,zs0] = zscore(yROT(1:3:end));
netANGLE = trainNetwork(xROT(:,:,1,1:3:end),newY,layers,options);
%{
%%
reluLayer
    maxPooling2dLayer(2,'Stride',2)
    convolution2dLayer(3,3)
%}
%% setup for optimized network(s)
rateSchedule = optimizableVariable('LearnRateSchedule',{'piecewise','none'},'Type','categorical');
rateDropFactor = optimizableVariable('LearnRateDropFactor',[0,1],'Type','integer');
rateDropPeriod = optimizableVariable('LearnRateDropPeriod',[2,9],'Type','integer');
L2reg = optimizableVariable('L2Regularization',[.0001,.0008]);
initLearnRate = optimizableVariable('InitialLearnRate',[.0001,.1],'Transform','log');
layerSize = optimizableVariable('layerSize',[5,15],'Type','integer');
layerNumber = optimizableVariable('layerNumber',[5,15],'Type','integer');
momentum = optimizableVariable('Momentum',[0,1]);
Bpara = [rateSchedule,rateDropFactor,rateDropPeriod,L2reg,initLearnRate,layerSize,layerNumber,momentum];
%% setup for angle network
xData = rotI;
xData = reshape(xData,[size(xData,1) size(xData,2) 1 size(xData,3)]);
yData = yROT;
[yData,zu0,zs0] = zscore(yData);
downSample = 4;
xData = xData(:,:,:,1:downSample:end);
yData = yData(1:downSample:end);
[trainIDX,testIDX] = dividerand(numel(yData),.8,.2,0);
maxTIME = 1.5*60*60;
maxE = 350;
func = cFlow('mySimpleCnn');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
[BayesObject_Angle,trainedNetwork_Angle] = func(Bpara,xData(:,:,:,trainIDX),yData(trainIDX),xData(:,:,:,testIDX),yData(testIDX),maxTIME,'gpu',maxE);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%% local train on angle
aP = cFlowLoader(BayesObject_Angle);
para = aP.XAtMinEstimatedObjective;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make options for training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = trainingOptions('sgdm',...
    'LearnRateSchedule',char(para.LearnRateSchedule),...
    'InitialLearnRate',para.InitialLearnRate,...
    'LearnRateDropFactor',para.LearnRateDropFactor,...
    'LearnRateDropPeriod',para.LearnRateDropPeriod,...
    'MaxEpochs',maxE,...
    'L2Regularization',para.L2Regularization,...
    'Momentum',para.Momentum,...
    'Plots','training-progress',...
    'Verbose',true,...
    'ExecutionEnvironment','parallel');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layers = [ ...
    imageInputLayer([size(xData,1) size(xData,2) 1])
    convolution2dLayer([para.layerSize para.layerSize],para.layerNumber)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% train network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trainedNetwork_Angle = trainNetwork(xData,yData,layers,options);
%% setup for maj
xData = SAM1;
xData = reshape(xData,[size(xData,1) size(xData,2) 1 size(xData,3)]);
yData = yS1;
[yData,zu1,zs1] = zscore(yData);
downSample = 4;
xData = xData(:,:,:,1:downSample:end);
yData = yData(1:downSample:end);
[trainIDX,testIDX] = dividerand(numel(yData),.8,.2,0);
maxTIME = 1.5*60*60;
maxE = 350;
func = cFlow('mySimpleCnn');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
[BayesObject_Major,trainedNetwork_Major] = func(Bpara,xData(:,:,:,trainIDX),yData(trainIDX),xData(:,:,:,testIDX),yData(testIDX),maxTIME,'gpu',maxE);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%% setup for minor axis
xData = SAM2;
xData = reshape(xData,[size(xData,1) size(xData,2) 1 size(xData,3)]);
yData = yS2;
[yData,zu2,zs2] = zscore(yData);
downSample = 4;
xData = xData(:,:,:,1:downSample:end);
yData = yData(1:downSample:end);
[trainIDX,testIDX] = dividerand(numel(yData),.8,.2,0);
maxTIME = 1.5*60*60;
maxE = 350;
func = cFlow('mySimpleCnn');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
[BayesObject_Minor,trainedNetwork_Minor] = func(Bpara,xData(:,:,:,trainIDX),yData(trainIDX),xData(:,:,:,testIDX),yData(testIDX),maxTIME,'gpu',maxE);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%% major
rr = SAM1;
layers = [ ...
    imageInputLayer([size(rr,1) size(rr,2) 1])
    convolution2dLayer([7 7],5)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','InitialLearnRate',0.01, 'MaxEpochs',150,'Plots','training-progress','ExecutionEnvironment','parallel');
xROT = reshape(rr,[size(rr,1) size(rr,2) 1 size(rr,3)]);
[newY,zu1,zs1] = zscore(yS1');
netMAJOR = trainNetwork(xROT,newY,layers,options);
% minor
rr = SAM2;
layers = [ ...
    imageInputLayer([size(rr,1) size(rr,2) 1])
    convolution2dLayer(11,11)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','InitialLearnRate',0.001, 'MaxEpochs',150,'Plots','training-progress','ExecutionEnvironment','parallel');
xROT = reshape(rr,[size(rr,1) size(rr,2) 1 size(rr,3)]);
[newY,zu2,zs2] = zscore(yS2');
netMINOR = trainNetwork(xROT,newY,layers,options);

%% display angle
close all
rr = rotI;
xROT = reshape(rr,[size(rr,1) size(rr,2) 1 size(rr,3)]);
close all
eP = 2901;
yP = trainedNetwork_Angle.predict(xROT(:,:,1,eP:(eP+99)));
plot(yP*zs0+zu0,'r')
hold on
plot(yROT(eP:(eP+99)),'k');
%% display MAJOR
rr = SAM1;
xROT = reshape(rr,[size(rr,1) size(rr,2) 1 size(rr,3)]);
close all
eP = 2900;
yP = netMAJOR.predict(xROT(:,:,1,eP:(eP+100)));
plot(yP*zs1+zu1,'r')
hold on
plot(yS1(eP:(eP+100)),'k');
%% display MINOR
rr = SAM2;
xROT = reshape(rr,[size(rr,1) size(rr,2) 1 size(rr,3)]);
close all
eP = 4;
yP = netMINOR.predict(xROT(:,:,1,eP*(1:100)));
plot(yP*zs2+zu2,'r')
hold on
plot(yS2(eP*(1:100)),'k');
%% whole or ALL
for e = 1:5
    test = imread(iFileList{e});
    angleFunc = @(I,P)getAngle(I,P,trainedNetwork_Angle,[zu0 zs0],netMAJOR,[zu1 zs1],netMINOR,[zu2 zs2],25,[]);
    plFunc = @(I)generateImageDomain(I,25);
    result{e} = applyFuncToLocation(test,angleFunc,plFunc,0,1);
end
%% say what
e=1
pickN = pick + 25;

test = imread(iFileList{e});
IDX = sub2ind(size(test),pickN(:,1),pickN(:,2));
angleFunc = @(I,P)getAngle(I,P,trainedNetwork_Angle,[zu0 zs0],netMAJOR,[zu1 zs1],netMINOR,[zu2 zs2],25,[]);
plFunc = @(I)IDX;
result{e} = applyFuncToLocation(test,angleFunc,plFunc,0,1);
%% decompose the results
toD = [];
for e = 1:numel(result)
    toD = [toD;squeeze(result{e})];
    result{e} = [];
    e
end
%%
%toD(:,end-2:end) = [];
[sU,sE] = PCA_FIT_FULLws(toD,9);
sC = PCA_REPROJ(toD,sE,sU);
%%
P = prod(size(test));
str = 1;
stp = str + P -1;
img = sC(str:stp,:);
img = reshape(img,[size(test) size(img,2)]);
vw = 1:3;
close all
for k = 1:size(img,3)
    img(:,:,k) = bindVec(img(:,:,k));
end
imshow(img(:,:,vw),[]);
waitforbuttonpress
close all
for k = 1:size(img,3)
    vw = k:(k+2);
    imshow(img(:,:,vw),[]);
    waitforbuttonpress
end
%%

%% view low dim from result
e = 5;
close all
para = squeeze(result{e}(:,:,end-2:end));
para = reshape(para,[size(test,1) size(test,2) 3]);
imshow(para(:,:,1),[])
%%
%result = reshape(squeeze(result),[size(test,1)-50 size(test,2)-50 3]);
for e = 1:4
    test(1:25,:) = [];
    test = imrotate(test,90);
end
%% view angle over point training data
close all
for e = 1:30%7%2%:numel(maizeStomataCenter_row)
    angleFunc = @(I,P)getAngle(I,P,trainedNetwork_Angle,[zu0 zs0],netMAJOR,[zu1 zs1],netMINOR,[zu2 zs2],25);
   
    tmpI = imread(iFileList{e});
    tmpP = [maizeStomataCenter_row{e},maizeStomataCenter_column{e}];
    %tmpP = new{e}(30,:) + [2 1];
    ridx = find(any(tmpP(:,1) < R | tmpP(:,2) < R | tmpP(:,1) > (size(tmpI,1) - R) | tmpP(:,2) > (size(tmpI,2) - R),2));
    tmpP(ridx,:) = [];
    tmpPI = sub2ind(size(tmpI),tmpP(:,2),tmpP(:,1));
    plFunc = @(N)tmpPI;
    imshow(tmpI,[]);
    drawnow
    hold on
    plot(tmpP(:,1),tmpP(:,2),'r*')
    probeD = applyFuncToLocation(tmpI,angleFunc,plFunc,0,1);
    TH = linspace(-pi,pi,100);
    dX10 = [cos(TH)' sin(TH)'];
    

    
    for p = 1:size(probeD,1)
        para = probeD(p,end-2:end);
        para = [para 0 0];
        para(1) = -para(1);
        [T] = generateDomainTransformation(para);
        T;
        c = tmpP(p,1);
        r = tmpP(p,2);
        T(:,3) = T(:,3) + [c;r];
        [dX1] = transformDomain(dX10,T);


        plot(dX1(:,1),dX1(:,2),'r');
        text(c-40,r-40,num2str(p),'BackgroundColor','w')
    end
    hold off
    drawnow
    pause(1)
title(num2str(e))
waitforbuttonpress
end

%% blah
close all
viewC = wC(:,:,1:3);
for e = 1:size(viewC,3)
    viewC(:,:,e) = bindVec(viewC(:,:,e));
end
imshow(viewC,[])
%%
close all
kW = reshape(wC,[size(wC,1)*size(wC,2) size(wC,3)]);
kidx = kmeans(kW,10);
RGB = reshape(kidx,[size(wC,1) size(wC,2)]);
RGB = label2rgb(RGB);
imshow(RGB,[]);

%% select process
close all
for e = 1
    test = imread(iFileList{e});
    test = imcrop(test);
    angleFunc = @(I,P)getAngle(I,P,trainedNetwork_Angle,[zu0 zs0],netMAJOR,[zu1 zs1],netMINOR,[zu2 zs2],25);
    plFunc = @(I)generateImageDomain(I,25);
    TESTER = applyFuncToLocation(test,angleFunc,plFunc,0,1);
end
for e = 1:4
    test(1:25,:) = [];
    test = imrotate(test,90);
end
TESTER = squeeze(TESTER(:,:,end-2:end));
TESTER = reshape(TESTER,[size(test,1) size(test,2) 3]);
%% final search
e=5;
cach = [];
cach2 = [];
RAW = [];
R = 25+4;
midx = [];

for e = 1:3

    testR = imread(iFileList{e});
    test = imread(iFileList{e});
    for r = 1:4
        test(1:25,:) = [];
        test = imrotate(test,90);
    end
    TESTER = squeeze(result{e}(:,end-3:end-1));
    TESTER = reshape(TESTER,[size(test) size(TESTER,2)]);

    tmpP = [maizeStomataCenter_row{e},maizeStomataCenter_column{e}];
    %tmpP = new{e}(30,:) + [2 1];
    ridx = find(any(tmpP(:,1) < R | tmpP(:,2) < R | tmpP(:,1) > (size(test,1) - R) | tmpP(:,2) > (size(test,2) - R),2));
    tmpP(ridx,:) = [];
    
    tmpP = tmpP - 25;
    tmpPI = sub2ind(size(test),tmpP(:,2),tmpP(:,1));
    plFunc = @(N)tmpPI;
    use = 'No';
    for p = 1:size(tmpP,1)
    %for p = 1:22
        dn = 'No';
        %while ~strcmp(dn,'Yes')
            BOX = [tmpP(p,:) 50 50];
            subI = imcrop(testR,BOX);
            subI = interp2(subI,2);
            %[se(1) se(2) V] = impixel(subI);
            
            se2 = round(se * .25-25);

            if size(cach,1) >= 2
                se2 = [0 0];
                dl = -4:4;
                Lc = 1;
                SAMP = [];
                for l1 = 1:numel(dl)
                    for l2 = 1:numel(dl)
                        para = squeeze(TESTER(se2(1)+tmpP(p,2)+dl(l1),se2(2)+tmpP(p,1)+dl(l2),:));
                        para = [para;0;0];
                        TH = linspace(-pi,pi,100);
                        dX10 = [cos(TH)' sin(TH)'];
                        [T] = generateDomainTransformation(para);
                        T(:,3) = T(:,3) + [tmpP(p,2)+dl(l1);tmpP(p,1)+dl(l2)];
                        IDX = sub2ind(size(test),se2(1)+tmpP(p,2)+dl(l1),se2(2)+tmpP(p,1)+dl(l2));
                        tmpSAMP = squeeze(result{e}(IDX,:,1:end-4));
                        tmpSAMP = zscore(tmpSAMP(:),1,1);
                        tmpSAMP = reshape(tmpSAMP,[50 50]);
                        tmpSAMP = cat(3,tmpSAMP,flipdim(tmpSAMP,2),flipdim(tmpSAMP,1),flipdim(flipdim(tmpSAMP,1),2));
                        tmpSAMP2 = std(tmpSAMP,1,3);
                        tmpSAMP = mean(tmpSAMP,3);
                        SAMP(Lc,:) = tmpSAMP(:)';
                        SAMP2(Lc,:) = tmpSAMP2(:)';
                        Lc = Lc + 1;
                    end
                end
                PROB = [];
                PROB2 = [];
                J1 = mean(cach,1);
                J2 = std(cach,1,1);
                K1 = mean(cach2,1);
                K2 = std(cach2,1,1);

                for i = 1:size(SAMP,1)
                    for j = 1:size(SAMP,2)
                        PROB(i,j) = normpdf(SAMP(i,j),J1(j),J2(j));
                        PROB2(i,j) = normpdf(SAMP2(i,j),K1(j),K2(j));
                    end
                    i
                end
                DELTA = sum(-log(PROB),2);
                DELTA2 = sum(-log(PROB2),2);
                DELTA = DELTA+DELTA2;
                %DELTA = bsxfun(@minus,SAMP,mean(cach,1));
                %DELTA = sum(DELTA.*DELTA,2);
                [~,midx] = min(DELTA);
            end
          

            se2 = [0 0];
            dl = -4:4;
            Lc = 1;
            if ~isempty(midx)
                [l2,l1] = ind2sub([numel(dl) numel(dl)],midx);
                para = squeeze(TESTER(se2(1)+tmpP(p,2)+dl(l1),se2(2)+tmpP(p,1)+dl(l2),:));
                para = [para;(25+se2(2)+dl(l2))*4;(25+se2(1)+dl(l1))*4];
                para2 = para;
                para2(2:3) = para2(2:3)*4;
                [T] = generateDomainTransformation(para2);
                [dX1] = transformDomain(dX10,T);
                imshow(subI,[]);
                hold on
                plot(dX1(:,1),dX1(:,2),'b');
                %use = questdlg('use?');
                %use = 'No';
            end
            if ~strcmp(use,'Yes')
                for l1 = 1:numel(dl)
                    for l2 = 1:numel(dl)
                        para = squeeze(TESTER(se2(1)+tmpP(p,2)+dl(l1),se2(2)+tmpP(p,1)+dl(l2),:));
                        para = [para;(25+se2(2)+dl(l2))*4;(25+se2(1)+dl(l1))*4];
                        para2 = para;
                        para2(2:3) = para2(2:3)*4;
                        TH = linspace(-pi,pi,100);
                        dX10 = [cos(TH)' sin(TH)'];
                        [T] = generateDomainTransformation(para2);
                        T
                        %para = [para;fliplr(se)'];
                        %T(:,3) = T(:,3) + [tmpP(p,1);tmpP(p,2)];
                        [dX1] = transformDomain(dX10,T);
                        imshow(subI,[]);
                        drawnow
                        hold on
                        %plot(tmpP(:,1),tmpP(:,2),'r*');
                        if Lc == midx
                            plot(dX1(:,1),dX1(:,2),'g');
                        else
                            plot(dX1(:,1),dX1(:,2),'r');
                        end


                        title(num2str(Lc));
                        hold off
                        Lc = Lc + 1;
                        drawnow
                        %waitforbuttonpress
                    end

                end
                %sel = inputdlg('Which?');
                %sel = str2num(sel{1});
            else
                sel = midx;
            end
            sel = SEL{e,p};
            %SEL{e,p} = sel;

            [l2 l1] = ind2sub([numel(dl) numel(dl)],sel);
            para = squeeze(TESTER(se2(1)+tmpP(p,2)+dl(l1),se2(2)+tmpP(p,1)+dl(l2),:));
            para = [para;(25+se2(2));(25+se2(1))];
            [T] = generateDomainTransformation(para);
            T(:,3) = T(:,3) + [tmpP(p,1)+dl(l2);tmpP(p,2)+dl(l1)];
            IDX = sub2ind(size(test),se2(1)+tmpP(p,2)+dl(l1),se2(2)+tmpP(p,1)+dl(l2));
            [dX1] = transformDomain(dX10,T);
            TT = imread(iFileList{e});
            close all
            imshow(TT,[]);
            drawnow
            hold on
            %plot(tmpP(:,1),tmpP(:,2),'r*');
            plot(dX1(:,1),dX1(:,2),'r');
            plot(se2(2)+tmpP(p,1)+dl(l2)+25,se2(1)+tmpP(p,2)+dl(l1)+25,'g.')
            drawnow
            tmpSAMP = squeeze(result{e}(IDX,:,1:end-4))';
            tmpSAMP = zscore(tmpSAMP(:),1,1);
            RAW = [RAW;tmpSAMP(:)'];
            tmpSAMP = reshape(tmpSAMP,[50 50]);
            tmpSAMP = cat(3,tmpSAMP,flipdim(tmpSAMP,2),flipdim(tmpSAMP,1),flipdim(flipdim(tmpSAMP,1),2));
            tmpSAMP2 = std(tmpSAMP,1,3);
            tmpSAMP = mean(tmpSAMP,3);
            cach = [cach;tmpSAMP(:)'];
            cach2 = [cach2;tmpSAMP2(:)'];
            %waitforbuttonpress
            %{
            angleFunc = @(I,P)getAngle(I,P,trainedNetwork_Angle,[zu0 zs0],netMAJOR,[zu1 zs1],netMINOR,[zu2 zs2],25);
            plFunc = @(I)IDX;
            what = applyFuncToLocation(test,angleFunc,plFunc,0,1);
            %}
            %dn = questdlg('done?');
        %end
        
    end
end
%% grade again
close all
[g1,g2] = ndgrid(linspace(-1.5,1.5,50),linspace(-1.5,1.5,50));
G = (g1.^2 + g2.^2).^.5;
imshow(G < 1);
gidxI = find(G<1);
gidxE = find(G>1);
ecach = [];
icach = [];
for e = 1:size(RAW,1)
    [d] = patchFlipOp(reshape(RAW(e,:),[50 50]));
    d = d';
    d = RAW(e,:);
    icach = [icach;d(:,gidxI)];
    ecach = [ecach;d(:,gidxE)];
%{
    icach(e,:) = RAW(e,gidxI);
    ecach(e,:) = RAW(e,gidxE);
    %}
    

    %icach2(e,:) = cach2(e,gidxI);
    %ecach2(e,:) = cach2(e,gidxE);
end

[icU,icE,icL]= PCA_FIT_FULLws(icach,3);
[ecU,ecE,ecL]= PCA_FIT_FULLws(ecach,3);
[icC]= PCA_REPROJ(icach,icE,icU);
[ecC]= PCA_REPROJ(ecach,ecE,ecU);
figure;
plot3(icC(:,1),icC(:,2),icC(:,3),'.')
figure;
plot3(ecC(:,1),ecC(:,2),ecC(:,3),'.')
figure;
[iSim]= PCA_BKPROJ(icC,icE,icU);
[eSim]= PCA_BKPROJ(ecC,ecE,ecU);
iErr = sum((iSim - icach).^2,2).^.5;
eErr = sum((eSim - ecach).^2,2).^.5;
iPDF_u = [zeros(1,size(icE,2)) mean(iErr)];
iPDF_c = diag([diag(icL);std(iErr)^2]);
ePDF_u = [zeros(1,size(ecE,2)) mean(eErr)];
ePDF_c = diag([diag(ecL);std(eErr)^2]);
maskI = zeros(size(g1));
maskI(gidxI) = icU;
imshow(maskI,[]);
maskE = zeros(size(g1));
maskE(gidxE) = ecU;
imshow(maskE,[]);
maskT = zeros(size(maskI));
maskT = maskI + maskE;
imshow(maskT,[]);
maskW = reshape(mean(cach,1),[50 50]);
imshow(maskW,[]);
MM = [];
for d = 1:size(icE,2)
    tM = [];
    sweep = linspace(-icL(d,d),icL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(icE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,icE,icU);
        tmpM = zeros(size(maskI));
        tmpM(gidxI) = tmpD;
        tM = [tM tmpM];
    end
    MM = [MM;tM];
end
imshow(MM,[]);
figure;
MM = [];
for d = 1:size(ecE,2)
    tM = [];
    sweep = linspace(-ecL(d,d),ecL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(ecE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,ecE,ecU);
        tmpM = zeros(size(maskI));
        tmpM(gidxE) = tmpD;
        tM = [tM tmpM];
    end
    MM = [MM;tM];
end
imshow(MM,[]);
% regress
rIE = icC\ecC;
rEI = ecC\icC;
ksI = [icC iErr];
ksE = [ecC eErr];
% ksdensity plots
for e = 1:size(ksI,2)
    figure
    [ipdff(e,:),ixi(e,:)] = ksdensity(ksI(:,e),'Bandwidth',.3*std(ksI(:,e)));
    ipdff(e,:) = ipdff(e,:)/sum(ipdff(e,:));
    plot(ixi(e,:),ipdff(e,:));
    title(['PC:' num2str(e) '<--internal']);
end
for e = 1:size(ksE,2)
    figure
    [epdff(e,:),exi(e,:)] = ksdensity(ksE(:,e),'Bandwidth',.3*std(ksE(:,e)));
    epdff(e,:) = epdff(e,:)/sum(epdff(e,:));
    plot(exi(e,:),epdff(e,:));
    title(['PC:' num2str(e) '<--external']);
end
% store PDF data
myiPDF.X = ixi;
myiPDF.F = ipdff;
myePDF.X = exi;
myePDF.F = epdff;
% regression plots
toScatterE = (icC*rIE);
for r = 1:size(toScatterE,2)
    figure;
    plot(ecC(:,r),toScatterE(:,r),'k.')
    hold on
    plot(linspace(min(ecC(:,r)),max(ecC(:,r)),10),linspace(min(ecC(:,r)),max(ecC(:,r)),10),'r');
end
% regression plots
toScatterI = (ecC*rEI);
for r = 1:size(toScatterI,2)
    figure;
    plot(icC(:,r),toScatterI(:,r),'k.')
    hold on
    plot(linspace(min(icC(:,r)),max(icC(:,r)),10),linspace(min(icC(:,r)),max(icC(:,r)),10),'r');
end
figure;
MM = [];
for d = 1:size(icE,2)
    tM = [];
    sweep = linspace(-icL(d,d),icL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(icE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,icE,icU);
        tmpM = zeros(size(maskI));
        tmpM(gidxI) = tmpD;
    
        tmpE = (tmpC*rIE);
        tmpE = PCA_BKPROJ(tmpE,ecE,ecU);
        tmpME = zeros(size(maskI));
        tmpME(gidxE) = tmpE;
        tmpM = tmpM + tmpME;

        tM = [tM tmpM];
    end
    MM = [MM;tM];
end

imshow(MM,[]);
title('internal->external');
drawnow
figure;
MM = [];
for d = 1:size(ecE,2)
    tM = [];
    sweep = linspace(-ecL(d,d),ecL(d,d),10);
    for s = 1:numel(sweep)
        tmpC = zeros(1,size(ecE,2));
        tmpC(d) = sweep(s);
        tmpD = PCA_BKPROJ(tmpC,ecE,ecU);
        tmpM = zeros(size(maskI));
        tmpM(gidxE) = tmpD;
    
        tmpI = (tmpC*rEI);
        tmpI = PCA_BKPROJ(tmpI,icE,icU);
        tmpMI = zeros(size(maskI));
        tmpMI(gidxI) = tmpI;
        tmpM = tmpM + tmpMI;

        tM = [tM tmpM];
    end
    MM = [MM;tM];
end

imshow(MM,[]);
title('external->internal');
%% fun grade
e = 1;
data = Gresult{e}(:,1:end-4);
[G] = issueGrades(data,gidxI,icE,icU,iPDF_u,iPDF_c,gidxE,ecE,ecU,ePDF_u,ePDF_c,rEI,rIE,myiPDF,myePDF);
%%
G = [G prod(G,2)];
%% peak for fun
e = 1;
tmpP = [maizeStomataCenter_row{e},maizeStomataCenter_column{e}];
ridx = find(any(tmpP(:,1) < R | tmpP(:,2) < R | tmpP(:,1) > (size(test,1) - R) | tmpP(:,2) > (size(test,2) - R),2));
tmpP(ridx,:) = [];

p=20;
close all
IDX = sub2ind(size(test),se2(1)+tmpP(p,2)+dl(l1)-25,se2(2)+tmpP(p,1)+dl(l2)-25);
%pull = result{e}(IDX,1:end-4);
%figure;
%pull = zscore(pull(:),1,1);
%pull = reshape(pull,[50 50]);
imshow(pull,[]);
[l2 l1] = ind2sub([numel(dl) numel(dl)],SEL{1,p});
para = squeeze(TESTER(se2(1)+tmpP(p,2)+dl(l1)-25,se2(2)+tmpP(p,1)+dl(l2)-25,:));
para = [para;0;0];
[T] = generateDomainTransformation(para);
T(:,3) = T(:,3) + [tmpP(p,1)+dl(l2);tmpP(p,2)+dl(l1)];
TH = linspace(-pi,pi,100);
dX10 = [cos(TH)' sin(TH)'];
[dX1] = transformDomain(dX10,T);
imshow(test,[]);
hold on
plot(dX1(:,1)-25,dX1(:,2)-25,'r')
plot(tmpP(p,1)-25+dl(l2),tmpP(p,2)-25+dl(l1),'r.');

%%
NAMES = {'internal' 'external' 'external-->internal' 'internal-->external' 'iter-internal' 'iter-external'};
e=1;
test = imread(iFileList{e});
for r = 1:4
    test(1:25,:) = [];
    test = imrotate(test,90);
end
%close all

%[hh1,hh2,~] = impixel(test);
hh2 = pick(2,1);
hh1 = pick(2,2);
TESTER = squeeze(Gresult{e}(:,end-3:end-1));
TESTER = reshape(TESTER,[size(test) size(TESTER,2)]);


for g = 1:size(G,2)
    tG = reshape(G(:,g),size(test));
    hope = -log(tG);
    msk = imerode(hope,strel('disk',15,0))==hope;
    msklarge = bwareaopen(msk,2);
    msk = msk == 1 & msklarge == 0;
    [hx1,hy1] = find(msk);
    figure;
    imshow(test,[]);
    hold on
    plot(hy1,hx1,'.');

    for p = 1:numel(hx1)
        para = squeeze(TESTER(hx1(p),hy1(p),:));
        para = [para;0;0];
        [T] = generateDomainTransformation(para);
        T(:,3) = T(:,3) + [hy1(p);hx1(p)];
        TH = linspace(-pi,pi,100);
        dX10 = [cos(TH)' sin(TH)'];
        [dX1] = transformDomain(dX10,T);
        plot(dX1(:,1),dX1(:,2),'r');
    end


    para = squeeze(TESTER(hh2,hh1,:));
    para = [para;0;0];
    [T] = generateDomainTransformation(para);
    T(:,3) = T(:,3) + [hh1;hh2];
    TH = linspace(-pi,pi,100);
    dX10 = [cos(TH)' sin(TH)'];
    [dX1] = transformDomain(dX10,T);
    plot(dX1(:,1),dX1(:,2),'g');
    plot(hh1,hh2,'g.');


    figure;
    imshow(tG,[]);
    waitforbuttonpress
    drawnow


    try
        title(NAMES{g})
    catch
        title('COMP')
    end
end
%%
eG = reshape(G(:,end),size(test));
tG = reshape(prod(G(:,end-1:end),2),size(test));
figure;
imshow(cat(3,bindVec(iG),bindVec(eG),bindVec(tG)),[]);
%%
pick(1,:) = [385 230];
pick(2,:) = [388 226];
IDX = sub2ind(size(test),pick(:,1),pick(:,2));
G(IDX,end)
TESTER(pick(1),pick(2),:)
[iG] = issueGrades(data(IDX,:),gidxI,icE,icU,iPDF_u,iPDF_c,gidxE,ecE,ecU,ePDF_u,ePDF_c,rEI,rIE,myiPDF,myePDF);
%%


hope = -log(eG);
[hx2,hy2] = find(imerode(hope,strel('disk',21,0))==hope);

hope = -log(tG);
[hx3,hy3] = find(imerode(hope,strel('disk',21,0))==hope);


figure;
pause(1);
imshow(test,[]);
hold on
plot(hy1,hx1,'r.')
plot(hy2,hx2,'g.')
plot(hy3,hx3,'b.')
%[hh1,hh2,~] = impixel();
[~,Q] = min(sum((hy1-hh1).^2 + (hx1-hh2).^2,2));
plot(hy1(Q),hx1(Q),'ro');
%IDX = sub2ind(size(test),hx1(Q),hy1(Q));
iG = issueGrades(data(IDX,:),gidxI,icE,icU,iPDF_u,iPDF_c,gidxE,ecE,ecU,ePDF_u,ePDF_c,rEI,rIE);
%% whole grade
[cU,cE,cL]= PCA_FIT_FULLws(cach,3);
junk = PCA_REPROJ(cach,cE,cU);
simJ = PCA_BKPROJ(junk,cE,cU);
errJ = sum((cach-simJ).^2,2).^.5;
mean(errJ)
funU = [zeros(1,size(junk,2)) mean(errJ)];
funC = diag([diag(cL);std(errJ).^2]);



[cU2,cE2,cL2]= PCA_FIT_FULLws(cach2,3);
junk2 = PCA_REPROJ(cach2,cE2,cU2);
simJ2 = PCA_BKPROJ(junk2,cE2,cU2);
errJ2 = sum((cach2-simJ2).^2,2).^.5;
funU2 = [zeros(1,size(junk2,2)) mean(errJ2)];
funC2 = diag([diag(cL2);std(errJ2).^2]);

grader.E = cE;
grader.U = cU;
grader.dU = funU;
grader.dC = funC;

grader.E2 = cE2;
grader.U2 = cU2;
grader.dU2 = funU2;
grader.dC2 = funC2;
%% whole or ALL
%gU = mean(cach,1);
%gS = std(cach,1,1);
for e = 1
    test = imread(iFileList{e});
    angleFunc = @(I,P)getAngle(I,P,trainedNetwork_Angle,[zu0 zs0],netMAJOR,[zu1 zs1],netMINOR,[zu2 zs2],25,grader);
    plFunc = @(I)generateImageDomain(I,25);
    Gresult{e} = applyFuncToLocation(test,angleFunc,plFunc,0,1);
end
%% hope
close all
e=1;
test = imread(iFileList{e});
for r = 1:4
    test(1:25,:) = [];
    test = imrotate(test,90);
end
hope = squeeze(Gresult{e});
hope = reshape(hope(:,end),[size(test)]);
[hx,hy] = find(imerode(hope,strel('disk',21),0)==hope);
imshow(test,[]);
hold on
plot(hy,hx,'r*');
%% select view
e=5;
test = imread(iFileList{e});
for r = 1:4
    test(1:25,:) = [];
    test = imrotate(test,90);
end

dl = -3:3;
TESTER = para;
n = 5;
str = (e-1)*n + 1;
stp = str + prod(size(test))-1;

subD = toD(str:stp,:);
TESTER = squeeze(result{e}(:,end-2:end));
TESTER = reshape(TESTER,[size(test) size(TESTER,3)]);
%%
e=1
TESTER = squeeze(Gresult{e}(:,end-3:end-1));
TESTER = reshape(TESTER,[size(test) size(TESTER,2)]);
%%
[co ro V] = impixel(test);
delta = sum((co-hy).^2 + (ro-hx).^2,2);
[~,mmidx ] = min(delta);
sresult = imfilter(TESTER,fspecial('gaussian',[21 21],2),'replicate');
dl = -4:4
for l1 = 1:numel(dl)
    for l2 = 1:numel(dl)
        r = ro + dl(l1);
        c = co + dl(l2);
        para = squeeze(TESTER(r,c,:));
        para(4:5) = 0;
        para(1) = -para(1);
        %para(1) = 0;
        %para(1) = -pi/8;
        %para(2) = 20;
        TH = linspace(-pi,pi,100);
        
        dX10 = [cos(TH)' sin(TH)'];
        [T] = generateDomainTransformation(para);
        T
        T(:,3) = T(:,3) + [c;r];
        [dX1] = transformDomain(dX10,T);
        imshow(test,[]);
        hold on
        plot(dX1(:,1),dX1(:,2),'r')
        plot(dX1(1,1),dX1(1,2),'r*')
        plot(hy,hx,'.');


        para = squeeze(TESTER(hx(mmidx),hy(mmidx),:));
        para(4:5) = 0;
        para(1) = -para(1);
        TH = linspace(-pi,pi,100);
        dX10 = [cos(TH)' sin(TH)'];
        [T] = generateDomainTransformation(para);
        T(:,3) = T(:,3) + [hy(mmidx);hx(mmidx)];
        [dX1] = transformDomain(dX10,T);
        hold on
        plot(dX1(:,1),dX1(:,2),'y')
        %{
        dX20 = [cos(TH)' sin(TH)'];
        [T] = generateDomainTransformation(para);
        T
        %T(:,3) = T(:,3) + [c;r];
        [dX2] = transformDomain(dX20,T);
        hold on
        plot(dX2(:,1),dX2(:,2),'g')
        
        plot(dX2(1,1),dX2(1,2),'g*')
        %}
        drawnow
        waitforbuttonpress
        hold off
    end
end

       
%%
sz = size(SAMP);
SSSSP = reshape(SAMP,[prod(sz(1:2)) sz(3)]);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL_T(SSSSP,1);
%UI = mean(SAMP,2);
%EI = zeros(size(UI));
%% blah - da - blah
x = {};
for e = 1:numel(maizeStomataCenter_row)
    blah{e} = [[maizeStomataCenter_row{e}],[maizeStomataCenter_column{e}]];
    iP = sub2ind(size(tmpI),blah{e}(:,2),blah{e}(:,1));
    tmpI = imread(iFileList{e});
    iPara = mean(storePARA,1);
    [t1 t2] = ndgrid(linspace(-2,2,200),linspace(-2,2,200));
    tD = [t2(:) t1(:)];
    %tmpI = imread(iFileList{e});
    r = [];
    options = optimoptions(@fmincon,'Display','iter','FiniteDifferenceType','central','FiniteDifferenceStepSize',[.5 .5 .5 .5 .5] ,'TypicalX',[1*pi/180 1 1 1 1]);
    parfor p = 1:numel(iP)
        minFunc = @(ip)phunnyML(I,iP(p),tD,ip,size(t1),sU,sE);
        r(p,:) = fmincon(minFunc,iPara,[],[],[],[],[-5*pi/180 5 5 -5 -5],[-5*pi/180 20 15 5 5],[],options);
    end
    x{e} = r;
    imshow(tmpI,[]);
    TH = linspace(-pi,pi,200);
    dx = [cos(TH)' sin(TH)' ones(numel(TH),1)];
    hold on
    for p = 1:numel(iP)
        [T] = generateDomainTransformation(x{e}(p,:));
        T(:,3) =  T(:,3) + blah{e}(p,:)';
        dY = (T*dx')';
        tP = mean(dY,1);
        plot(dY(:,1),dY(:,2),'r');
        text(tP(1),tP(2)-20,num2str(p),'BackgroundColor','w')
        plot(blah{e}(p,1),blah{e}(p,2),'r*')
    end
    drawnow
end
%%

%% color images from clicks
I = imread(iFileList{1});
border = 40;
close all
Y = [];
for e = 1:numel(maizeStomataCenter_row)
    tmp = zeros(size(I));
    tmpI = imread(iFileList{e});
    for p = 1:size(maizeStomataCenter_row{e},1)
         tmp(maizeStomataCenter_column{e}(p),maizeStomataCenter_row{e}(p)) = 1;
    end
    tmp = imdilate(tmp,strel('disk',9,0));
    %tmp = tmp;
    tmp2 = imdilate(tmp,strel('disk',4,0));
    tmp = tmp + tmp2;
    for rot = 1:4
        tmp(1:border,:) = [];
        tmpI(1:border,:) = [];
        tmp = imrotate(tmp,90);
        tmpI = imrotate(tmpI,90);
    end
    %tmp = padarray(tmp,[1 1],0,'pre');
    %tmpI = padarray(tmpI,[1 1],0,'pre');
    %size(tmp)
    
    %out = flattenMaskOverlay(tmpI,logical(tmp));
    %imshow(out,[]);
    Y = [Y;tmp(:)];
    
    
    
    drawnow
end
close all
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Color images - end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract images - start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% extract bug eye
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create sample disk
tmpI = imread(iFileList{1});
R = [0 40];
N = [(R(2)-R(1)) round(2*pi*(R(2)))];
[n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
[d1,d2] = ndgrid((R(2)+1):(size(tmpI,1)-R(2)),(R(2)+1):(size(tmpI,2)-R(2)));
indexPosition = sub2ind(size(tmpI),d1(:),d2(:));
Xd = n1.*cos(n2);
Yd = n1.*sin(n2);
Domain = [Xd(:) Yd(:)];
% create position index list function
plFunc = @(I)generateImageDomain(I,R(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplerFunction = @(I,P)myInterp2Sampler(I,P,Domain,N);
funcToApply = @(I,P)fftPatch(I,P,samplerFunction,[1:15]);
for e = 1:30
    fprintf(['start:' num2str(e) '\n']);tic
    tmpI = double(imread(iFileList{e}))/255;
    [fft{e}] = applyFuncToLocation(tmpI,funcToApply,plFunc,0,0);
    fprintf(['end:' num2str(e) ':' num2str(toc) '\n'])
end
%% stack data
imgTOT = 433*433;
imgTOT = (size(fft{1},1)^.5)^2;
RAD = size(fft{1},2);
TH = size(fft{1},3);
X = zeros([TH RAD 1 imgTOT*numel(fft) size(fft{1},4)]);
str = 1;
for e = 1:numel(fft)
    %for stackDim = 1:size(fft{1},4)
        fprintf(['start stacking:' num2str(e) ':' num2str(numel(fft)) '\n']);
        tmp = permute(fft{e},[3 2 1 4]);
        stp = str + size(tmp,3) - 1;
        X(:,:,1,str:stp,:) = tmp;
        str =  stp + 1;
        fprintf(['end stacking:' num2str(e) ':' num2str(numel(fft)) '\n']);
    %end
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
%% reduce dims along radius - PCA only the cliked points
sidx = find(double(cY)==3);
for d = 1:size(X,5)
    fX = squeeze(X(:,:,1,sidx,d));
    fX = reshape(fX,[size(fX,1) size(fX,2)*size(fX,3)]);
    numChoose = 10;
    %[qS qC qU qE qL qERR qLAM] = PCA_FIT_FULL_T(fX,size(fX,1));
    [newU{d},newE{d}] = PCA_FIT_FULL_Tws(fX,numChoose);


    fX = squeeze(X(:,:,1,:,d));
    fX = reshape(fX,[size(fX,1) size(fX,2)*size(fX,3)]);
    newC = PCA_REPROJ_T(fX,newE{d},newU{d});
    simD = PCA_BKPROJ_T(newC,newE{d},newU{d});
    errD = sum((simD - fX).^2,1).^.5;
    newC = reshape([newC;errD],[numChoose+1 size(X,2) size(X,4)]);
    outC{d} = newC;
end
%%
newC = [outC{1} ; outC{2}];
%% LSTN
SIG = cell(1,size(X,4));
for e = 1:size(X,4)
    %SIG{e} = X(:,:,1,e);
    %SIG{e} = newC(:,:,e);
    SIG{e} = newC(:,:,e);
    (e/size(X,4))
end
%% make categorical array
cY = categorical(Y);
%% mix proper ratio of ones and zeros
fidx0 = find(Y==0);
fidx1 = find(Y==1);
fidx2 = find(Y==2);

fidxM = [fidx1(randperm(round(numel(fidx1))));fidx2(randperm(round(numel(fidx2))));fidx0(randperm(round(numel(fidx1))))];
fidxM = fidxM(randperm(numel(fidxM)));
%% LSTN
options = trainingOptions('sgdm',...
            'InitialLearnRate',.003,...
            'MaxEpochs',1000,...
            'Shuffle','every-epoch',...
            'Verbose',true,...
            'ExecutionEnvironment','cpu',...
            'Plots','training-progress');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layers = [ ...
    sequenceInputLayer(size(SIG{1},1))
    lstmLayer(7,'OutputMode','last')
    fullyConnectedLayer(numel(unique(cY)))
    softmaxLayer
    classificationLayer];

%trainedNet_FM = trainNetwork(SIGGGY,YY,layers,options);

trainedNet_FM = trainNetwork(SIG(fidxM),cY(fidxM),layers,options);
%%
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplerFunction = @(I,P)myInterp2Sampler(I,P,Domain,N);
preprocessFunction = @(I,P)fftPatch(I,P,samplerFunction,[1:15]);
funcToApply = @(I,P)applyNetworkToPoint(I,P,preprocessFunction,trainedNet_FM,newE,newU);
for e = 1
    fprintf(['start:' num2str(e) '\n']);tic
    
    tmpI = double(imread(iFileList{e}))/255;
    [pre{e}] = applyFuncToLocation(tmpI,funcToApply,plFunc,0,1);
    
    
    fprintf(['end:' num2str(e) ':' num2str(toc) '\n'])
end
%% call stomata master
clear applyNetworkToPoint plFunc preprocessFunction
tmpI = imread(iFileList{1});
% define the radius to look over
R = [0 40];
% define the number of points
N = [(R(2)-R(1)) round(2*pi*(R(2)))];
% create the sample grid
[n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
% create the apply grid
[d1,d2] = ndgrid((R(2)+1):(size(tmpI,1)-R(2)),(R(2)+1):(size(tmpI,2)-R(2)));
% generate the index position(s) of the apply grid
indexPosition = sub2ind(size(tmpI),d1(:),d2(:));
% generate the sample patch
Xd = n1.*cos(n2);
Yd = n1.*sin(n2);
% store the domain
Domain = [Xd(:) Yd(:)];
% create position index list function
%plFunc = @(I)(1:numel(I));
%plFunc = @(I)indexPosition;
plFunc = @(I)generateImageDomain(I,R(2));
% generate the sampler function
samplerFunction = @(I,P)myInterp2Sampler(I,P,Domain,N);
% generate the preprocess function
preprocessFunction = @(I,P)fftPatch(I,P,samplerFunction,[1:15]);
% generate the network apply
networkApply = @(I,P)applyNetworkToPoint(I,P,preprocessFunction,trainedNet_FM2,newE,newU);

newSZ = size(d1);
%%
stomataMaster(iFileList{70},networkApply,plFunc,newSZ);
%%
stomataMaster(iFileList{60},trainedNet_FM2,R,Domain,N,1:15,newE,newU,newSZ,'./output/','');
%% publish stomataMaster
func = @(FN)stomataMaster(FN,networkApply,plFunc,newSZ,'./output/','');
pF = partialFunction(func,'maizeStomataNNapp');
pF.publish();
%%
SLAP = 750;
SLAP = 600;
TOTE = round(.7*numel(fidxM));
SLAP = 1:TOTE;
SLAP2 = TOTE:numel(fidxM);
func = cFlow('hyperPdeploy_emerge');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
% max function evaluations
maxE = 20;
% slow down training to get stable error
toSlow = 1;
% 
maxEval = 10000;
% max time for training hyper parameters
maxTime = 6*60*60;
% execution environment gpu
exeEnvironment = {'gpu',false};
% setup the variables to optimize
optimVars = [
    optimizableVariable('lstmStates',[2 10],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
% call the func and store the result in beta0 = b0
b_sto2 = func(optimVars,SIG(fidxM(SLAP)),...
                           cY(fidxM(SLAP)),...
                            SIG(fidxM(SLAP2)),...
                            cY(fidxM(SLAP2)),...
                            exeEnvironment,maxE,toSlow,maxE,maxTime,size(SIG{1},1),'last');

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%% LOCAL EMERGE
hPara = cFlowLoader(b_sto2);
%hPara.XAtMinObjective = hPara.XAtMinObjective*100;
[valErrorE3,consE3,trainedNet_FM2] = makeObjFcn_emerge(...
                                    hPara.XAtMinEstimatedObjective,...
                                    SIG(fidxM),...
                                    cY(fidxM),...
                                    '',...
                                    '',...
                                    0,...
                                    'training-progress',...
                                    2000,...
                                    'cpu',...
                                    size(SIG{1},1),...
                                    'last');
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
      
layers = [imageInputLayer([size(X,1) size(X,2) 1])
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


[BayesObject] = hyperPdeploy(optimVars,...
                            X(:,:,:,find(cTR.test(1))),...
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

            
