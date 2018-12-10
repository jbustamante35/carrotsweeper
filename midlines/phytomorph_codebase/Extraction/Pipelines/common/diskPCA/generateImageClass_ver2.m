function [boxNETS] = generateImageClass_ver2(rootPath,inMEM)
    %% notes:
    %% this pipeline was put together while playing nerf gun fighting with my son Jordan and his friend
    %% ben bremer.  i was working hard at coding and trying to structure the code but was having problems.
    %% during tis time, the kids were asking me to play nerf gun fighting.  while i was resistant at first, because I 
    %% "needed" to get some stuff done, i chosse to take a break and play a few short rounds.  while puping nerf bullets
    %% though the gun, i began to build up the idea that i could look at the data as a flow through a spray nozzle.  each data point
    %% can now be seen as a "point" in a 1-order tensor space.  if this view is useful, then this will be extended to process 
    %% n-tensor sprays.  the perspective is meant to process data-streams and the players are data sources, nozzles, nozzle-manifolds.
    %% results:
    %% the flow is a lazy style execution.  the objects created are often functions that are meant to process the data
    %% rather than directly coding objects that proces the data.  this is not to say that i do NOT need to code
    %% the functions. but rather the "cartrigaes" that are coded at loaded into the nozzle etc. thus giving them common interfaces.
    %% objects list:
    %% reduction nozzle: this is a nozzle that reduces the data spray from N to M dimensions.  it can be attached to any nozzle
    %% or data source if it is compatable.
    %% transduction nozzle: this nozzle emitts a new data type after the data is pulled through it.  the function that creates this
    %% takes a function that returns a function.
    %% the resulting 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scan base path
    cdir = dir(rootPath);
    cdir(1:2) = [];
    masterList = {};
    typeList = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the image lists
    for e = 1:numel(cdir)
        FilePath = [rootPath '/' cdir(e).name '/imageStacks/'];
        FileList = {};
        FileExt = {'tif','TIF','nef','NEF','PNG','png','jpeg','JPG','jpg','tiff'};
        FileList = gdig(FilePath,FileList,FileExt,1);
        masterList = cat(2,masterList,FileList);
        typeList = cat(1,typeList,e*ones(numel(FileList),1));
    end
    
    
    
    thumb = zeros([50 50 3 numel(masterList)]);
    parfor e = 1:numel(masterList)
        tmpI = imread(masterList{e});
        tmpI = imresize(tmpI,[50 50]);
        thumb(:,:,:,e) = tmpI;
        e
    end
    
    numClasses = numel(unique(typeList));
    layers = [imageInputLayer([50 50 3]);
          convolution2dLayer([21 21],20);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(numClasses);
          softmaxLayer();
          classificationLayer()];
      
     options = trainingOptions('sgdm','MaxEpochs',20,...
                'InitialLearnRate',0.0001);

     typeNet = trainNetwork(thumb,categorical(typeList),layers,options);
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the image-list data source 
    imageSource = dataSource(masterList,@(X)double(imread(X)),1);
    % create nozzle for dataSource
    thumb_Nozzle = dataNozzle(@(X,e0,e1)func_thumbNail(X,[50 50],1,false),imageSource,1);
    % get the reduction nozzle for the thumb_nozzle
    rThumb_nozzle = thumb_Nozzle.getReductionNozzle(3);
    % transform nozzle into transduction nozzle
    % self-note: comment the transduction nozzle
    tThumb_nozzle = rThumb_nozzle.getTransductionNozzle(@(X)getTrainedNB(X,typeList));
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each type directory - type examples might be carrots, maize, etc
    for t = 1:numel(cdir)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the image list for the t-th type folder
        for e = 1:numel(cdir)
            FilePath = [rootPath '/' cdir(t).name '/imageStacks/'];
            FileList = {};
            FileExt = {'tif','TIF','nef','NEF','PNG','png','jpeg','JPG','jpg','tiff'};
            FileList = gdig(FilePath,FileList,FileExt,1);
        end
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the number of crop boxes for this type
        cropboxRoot = [rootPath '/' cdir(t).name '/cropBoxes/'];
        cropboxPath = dir(cropboxRoot);
        cropboxPath(1:2) = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each crop box type
        for c = 1:numel(cropboxPath)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the image lists
            for e = 1:numel(cdir)
                FilePath = [rootPath '/' cdir(t).name '/cropBoxes/' cropboxPath(c).name filesep];
                mFileList = {};
                FileExt = {'tif','TIF','nef','NEF','PNG','png','jpeg','JPG','jpg'};
                mFileList = gdig(FilePath,mFileList,FileExt,1);
            end
            
            
            
            % percent to reduce the image size
            reducPer = .5;
            perToTrain = .75;
            numToTrain = perToTrain*round(numel(FileList));
            numToTest = numel(FileList) - numToTrain;
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read training image stack - X
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['start reading train images \n']);
            [imageVecsTrain] = readImagesFromList(FileList(1:numToTrain),reducPer);
            fprintf(['stop reading train images \n']);
            
            fprintf(['start reading test images \n']);
            [imageVecsTest] = readImagesFromList(FileList((numToTrain+1):numel(FileList)),reducPer);
            fprintf(['stop reading test images \n']);
            
            fprintf(['start reading train masks \n']);
            [maskVecsTrain] = readImagesFromList(mFileList(1:numToTrain),reducPer);
            fprintf(['stop reading train masks \n']);
            
            fprintf(['start reading test masks \n']);
            [maskVecsTest] = readImagesFromList(mFileList((numToTrain+1):numel(FileList)),reducPer);
            fprintf(['stop reading test masks \n']);
            
            
            windowSize = 11;
            [Xtrain,Ytrain,YtrainSequence,YtrainSequence2] = stackSlidingData(imageVecsTrain,maskVecsTrain,windowSize);
            [Xtest,Ytest,YtestSequence,YtestSequence2] = stackSlidingData(imageVecsTest,maskVecsTest,windowSize);
           
            
            [U,E] = getCompressCropBoxMethod(Xtrain,3);
            [Xtrain] = compressForCropBoxMethod(Xtrain,U,E);
            [Xtest] = compressForCropBoxMethod(Xtest,U,E);
            
            
            
            [XtrainCOL,YtrainCOL,numClassesC] = stackColumnData(imageVecsTrain,maskVecsTrain,true);
            [XtestCOL,YtestCOL,numClassesCTE] = stackColumnData(imageVecsTest,maskVecsTest,true);
            
            [XtrainROW,YtrainROW,numClassesR] = stackRowReSizeData(imageVecsTrain,maskVecsTrain,200);
            [XtestROW,YtestROW,numClassesRTE] = stackRowReSizeData(imageVecsTest,maskVecsTest,200);
            
            
            
            
            [Uc,Ec] = getCompressCropBoxMethod(XtrainCOL,3);
            [XtrainCOL] = compressForCropBoxMethod(XtrainCOL,Uc,Ec);
            [XtestCOL] = compressForCropBoxMethod(XtestCOL,Uc,Ec);
            
            
            
            
            [Ur,Er] = getCompressCropBoxMethod(XtrainROW,3);
            [XtrainROW] = compressForCropBoxMethod(XtrainROW,Ur,Er);
            [XtestROW] = compressForCropBoxMethod(XtestROW,Ur,Er);
            
            
            
            
            % augment with location inforamtion 
            for e = 1:numel(Xtrain)
                Xtrain{e} = [Xtrain{e};(1:size(Xtrain{e},2))*size(Xtrain{e},2)^-1];
            end
            
            for e = 1:numel(Xtrain)
                Xtrain{e} = [Xtrain{e};(1:size(Xtrain{e},2))*size(Xtrain{e},2)^-1];
            end
            
            
            % train column sweep
            inputSize = size(XtrainCOL{1},1);
            outputSize = 9;
            outputMode = 'sequence';
            numClasses = numClassesC;
            layers = [ ...
                sequenceInputLayer(inputSize)
                lstmLayer(outputSize,'OutputMode',outputMode)
                fullyConnectedLayer(numClasses)
                softmaxLayer
                classificationLayer];
            
            maxEpochs = 1000;
            options = trainingOptions('sgdm', ...
                'InitialLearnRate',0.001,...
                'MaxEpochs',maxEpochs);
            netCOL = trainNetwork(XtrainCOL,YtrainCOL',layers,options);
            
            
            try
                % train column sweep
                inputSize = size(XtrainROW{1},1);
                outputSize = 5;
                outputMode = 'sequence';
                numClasses = numClassesR;
                layers = [ ...
                    sequenceInputLayer(inputSize)
                    lstmLayer(outputSize,'OutputMode',outputMode)
                    fullyConnectedLayer(numClasses)
                    softmaxLayer
                    classificationLayer];

                maxEpochs = 200;
                options = trainingOptions('sgdm', ...
                    'InitialLearnRate',0.001,...
                    'MaxEpochs',maxEpochs);
                netROW = trainNetwork(XtrainROW,YtrainROW',layers,options);
            catch ME
                ME
            end
            
            boxNETS.class = typeNet;
            boxNETS.compress{t}{c}.columnCompressMean  = Uc;
            boxNETS.compress{t}{c}.columnCompressFrame  = Ec;
            boxNETS.compress{t}{c}.rowCompressMean  = Ur;
            boxNETS.compress{t}{c}.rowCompressFrame  = Er;
            boxNETS.drawers{t}{c}.columnNet = netCOL;
            boxNETS.drawers{t}{c}.rowNet = netROW;            
            
            
            %{
            %% remote GPU - train on GPU via condor and hyper parameters
           
            func = cFlow('hyperPdeploy_classify');
            func.setMCRversion('v930');
            func.setMemory('8000');
            func.setGPU(1);
            
            
            % max function evaluations
            maxE = 1000;
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
                optimizableVariable('InitialLearnRate',[1e-5 5e-2],'Transform','log')
                optimizableVariable('Momentum',[0.8 0.95])
                optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
            % call the func and store the result in beta0 = b0
            hyperParameters = func(optimVars,XtrainCOL,...
                                       YtrainCOL',...
                                        XtestCOL,...
                                        YtestCOL',...
                                        exeEnvironment,maxE,toSlow,maxE,maxTime,size(XtrainCOL{1},1),'sequence');

            auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
            auth = auth{1};
            func.submitDag(auth,50,50);
            
            
            %% perform local training with hyper parameters from condor
            hPara = cFlowLoader(hyperParameters);
            [valError,cons,trainedNet] = makeObjFcn_classify(...
                                                hPara.XAtMinObjective,...
                                                XtrainCOL,...
                                                YtrainCOL',...
                                                '',...
                                                '',...
                                                0,...
                                                'training-progress',...
                                                200,...
                                                'cpu',...
                                                size(XtrainCOL{1},1),...
                                                'sequence');
            
            
            %}
            
            
            
            %{
            inputSize = size(Xtrain{1},1);
            outputSize = 5;
            outputMode = 'last';
            numClasses = 2;
            layers = [ ...
                sequenceInputLayer(inputSize)
                lstmLayer(outputSize,'OutputMode',outputMode)
                fullyConnectedLayer(numClasses)
                softmaxLayer
                classificationLayer];
            
            maxEpochs = 5;
            options = trainingOptions('sgdm', ...
                'InitialLearnRate',0.1,...
                'MaxEpochs',maxEpochs);
            net = trainNetwork(Xtrain,Ytrain',layers,options);
            %}
            
            
            
            %{
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
            hyperParameters = func(optimVars,Xtrain,...
                                       Ytrain',...
                                        Xtest,...
                                        Ytest',...
                                        exeEnvironment,maxE,toSlow,maxE,maxTime,size(Xtrain{1},1),'last');

            auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
            auth = auth{1};
            func.submitDag(auth,50,50);
            
            
            %}
            
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TEST
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            toTest = 89;
            [testImageStack] = readImagesFromList(FileList(toTest),reducPer);
            
            
            [testImageStackCOL] = stackColumnData(testImageStack,maskVecsTrain);
            [testImageStackCOL] = compressForCropBoxMethod(testImageStackCOL,Uc,Ec);
            
            
            %{
            [testImageStack,~,~] = stackSlidingData(testImageStack,[],windowSize);
            [testImageStack] = compressForCropBoxMethod(testImageStack,U,E);
            for e = 1:numel(testImageStack)
                testImageStack{e} = [testImageStack{e};(1:size(testImageStack{e},2))*size(testImageStack{e},2)^-1];
            end
            %}
            
            
            [predCOL scoreCOL] = classify(netCOL,testImageStackCOL);
            pulseCOL = ~(double(mod(double(predCOL{1}),2))==1);
            testI = imread(FileList{toTest});
            pulseCOL_WHOLE = imresize(double(pulseCOL),[1 size(testI,2)],'nearest');
            wholeMSK = repmat(pulseCOL_WHOLE,[size(testI,1) 1]);
            reMSK = imresize(wholeMSK,[size(testImageStack,1) size(testImageStack,2)],'nearest');
            
            
            [testImageStackROW] = stackRowReSizeData(testImageStack,reMSK,200);
            [testImageStackROW] = compressForCropBoxMethod(testImageStackROW,Ur,Er);
            [predROW scoreROW] = classify(netROW,testImageStackROW);
            pulseROW = ~(double(mod(double(predCOL{1}),2))==1);
            pulseMSK = double(pulseROW)'*double(pulseCOL);
            
            pulseMSK = imresize(pulseMSK,[size(testI,1) size(testI,2)],'nearest');
            
            
            
            out = flattenMaskOverlay(testI,logical(pulseMSK==1));
            imshow(out,[])
            %{
            [YPred score] = classify(net,testImageStack);
            testI = imread(FileList{toTest});
            YPred = double(YPred)-1;
            pulse = imresize(double(YPred'),[1 size(testI,2)],'nearest');
            pulse = repmat(pulse,[size(testI,1) 1]);
            out = flattenMaskOverlay(testI,logical(pulse==1));
            %}
            
            
            
            
            
            
            
            
            %}
            
        end
    end
end
    
    
    
     
    %{
    
    
    vertStripMask.resetPtr();
    vertStrip.resetPtr();
    for e = 1:10
        I = vertStrip.next();
        
        M = vertStripMask().next();
        out = flattenMaskOverlay(I/280,logical(M));
        imshow(out,[]);
        drawnow
    end
    
    % spray test
    verticalANYnozzle.resetPtr();
    tester.resetPtr();
    cnt = 1;
    while tester.hasNext()
        gidx = tester.next();
        tidx = verticalANYnozzle.next();
        I = SUBimageSource.read(cnt);
        gidx = imresize(gidx,[1 size(I,2)]);
        tidx = imresize(tidx,[1 size(I,2)]);
        msk = repmat(gidx,[size(I,1) 1]);
        msk2 = repmat(tidx,[size(I,1) 1]);
        out=  flattenMaskOverlay(double(I)/255,logical(msk));
        out = flattenMaskOverlay(out,logical(msk2),.1,'b');
        imshow(out,[]);
        drawnow
        %{
        plot(gidx,'r');
        hold on
        plot(tidx,'b')
        drawnow
        hold off
        %}
        cnt = cnt + 1;
    end
    %}
    %{
    
    vertAny.resetPtr();
    tester2.resetPtr();
    
    cnt = 1;
    while tester2.hasNext()
        gidx = tester2.next();
        tidx = vertAny.next();
        I = vertStripImage.read(cnt);
        I = reshape(I,[size(I,1)/3 3 size(I,2)]);
        I = permute(I,[1 3 2]);
        gidx = imresize(gidx,[1 size(I,2)]);
        tidx = imresize(tidx,[1 size(I,2)]);
        msk = repmat(gidx,[size(I,1) 1]);
        msk2 = repmat(tidx,[size(I,1) 1]);
        I = permute(I,[2 1 3]);
        out=  flattenMaskOverlay((double(I)-min(I(:)))/350,logical(msk'));
        out = flattenMaskOverlay(out,logical(msk2'),.1,'b');
        imshow(out,[]);
        drawnow
        %{
        plot(gidx,'r');
        hold on
        plot(tidx,'b')
        drawnow
        hold off
        %}
        cnt = cnt + 1;
    end
    %}
    %{
    for e = 1:2
        filteredClassNozzle{e} = nozzleManifold({rSUBimageSource,verticalANYnozzle},@(X)filterSpray(X,e-1));
        meanGauge{e} = nozzleGauge(filteredClassNozzle{e},{@(X,v)particleCounter(X,v),@(X,v)accumulatorFunc(X,@(V)sum(V,2),v)});
        M{e} = meanGauge{e}.takeMetric();
        covGauge{e} = nozzleGauge(filteredClassNozzle{e},{@(X,v)particleCOV(X,M{e},v)});
        C{e} = covGauge{e}.takeMetric();
    end
    
    F = filteredClassNozzle{1}.collectFullSpray();
    
    
    
    
    [hmm] = makeChainRepeat(M{1}{2}/M{1}{1},C{1}{1},M{2}{2}/M{2}{1},C{2}{1},.999,3,nCOMP);
    
    func = @(X,e0,e1)hmm.Viterbi(X,ones(nCOMP,1),1);
    tester = dataNozzle(func,rSUBimageSource,1);
    tester.resetPtr();
    verticalANYnozzle.resetPtr();
    while tester.hasNext()
        gidx = tester.next();
        gidx = 1-mod(gidx,2);
        tidx = verticalANYnozzle.next();
        plot(gidx,'r')
        hold on
        plot(tidx,'b')
        drawnow
        hold off
    end
    
    
    
    
    meanGauge = nozzleGauge(thumb_Nozzle,{@(X,v)particleCounter(X,v),@(X,v)accumulatorFunc(X,@(V)sum(V,2),v)});
    M = meanGauge.takeMetric();
    covGauge = nozzleGauge(thumb_Nozzle,{@(X,v)particleCOV(X,M,v)});
    C = covGauge.takeMetric();
    
    
    
   
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create thumbnail-vec func
    toVecFunc = @(X,e)func_thumbNail(X,[50 50],1);
    % create select all func
    selIDFunc = @(X,e)X;
    % call PCA with source, thumbnail-vec, and select all
    [U,V] = diskPCA(I,toVecFunc,selIDFunc,5);
    % re-project to get comps
    [C] = diskBPROJ(I,U,V,toVecFunc,selIDFunc);
    
    
    
    % fit NB
    nb = fitcnb(C',typeList);
    
    
    
    
    
    
    
    
    
    
    UQ = unique(typeList);
    % create resize and depth stack func
    toVVecFunc = @(X,e)func_resizeDepthStack(X,.25,0);
    
    % number of components
    nCOMP = 4;
    % hardline return to node value
    HARDLINE_HOLD = .999;
    
    
    for u = 1:numel(UQ)
        
        % find the trials for type U
        fidx = find(UQ(u) == typeList);
        
        % call PCA on vertical stacked vectors
        [hU{u},hV{u},hD{u}] = diskPCA(I(fidx),toVVecFunc,selIDFunc,nCOMP);
        
       
       
        % load masks for this group
        tmpMsk = {};
        widthAverage = [];
        for e = 1:numel(fidx)
            fName = strrep(masterList{fidx(e)},'imageStacks','cropBoxes');
            fName = strrep(fName,'nef','tif');
            tmpMsk = imread(fName);
            tmpMsk = imresize(tmpMsk,per,'nearest');
            sig(e,:) = any(tmpMsk==1,1);
            R = regionprops(sig(e,:));
            widthAverage(e) = sum(sig(e,:));
            CNT(e) = numel(R);
        end
        widthAverage = sum(widthAverage)/(numel(widthAverage)*mean(CNT));
        
        
        
        
        
        
        toSelk0 = @(X,e)selectVec_ver0(X,sig,e,0);
        toSelk1 = @(X,e)selectVec_ver0(X,sig,e,1);
        getC = @(X,e)diskBPROJ({X},hU{u},hV{u},toVVecFunc,selIDFunc);
        [U0{u},C0{u}] = diskMEANandCOV(I(fidx),getC,toSelk0);
        [U1{u},C1{u}] = diskMEANandCOV(I(fidx),getC,toSelk1);
        [hmm{u}] = makeChainRepeat(U0{u},C0{u},U1{u},C1{u},HARDLINE_HOLD,mean(CNT),nCOMP);
        
       
        
        
        
        
        for e = 1:50
            gidx = hmm{u}.Viterbi(trainC{e},ones(nCOMP,1),1);
            gidx = 1-mod(gidx,2);
            close all
            plot(gidx);hold on;plot(sig(e,:),'r');
            drawnow
        end
        
        
        
        
        
        
        
    end
    %}
    
    

%{
    rootPath = '/mnt/snapper/nate/phunny/';
    T = generateImageClass_ver2(rootPath,1);
%}