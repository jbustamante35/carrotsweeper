function [T] = generateImageClass(rootPath,inMEM)
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
            
            
            
            % hard coded values for:
            % 1: percent resize 
            % 2: sort the data for markov-chain
            % 3: number of components for pls-regression
            % 4: transistion probs
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PER = .25;
            toSort = true;
            tPROB = .999;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % pseudo code
            % 1: construct the dataSource for images and masks
            % 2: construct data nozzle which resizes, sorts and depth
            % stacks the images
            % 3: constuct the data nozzle which resizes, sorts (false) and
            % depth stacks the masks
            % 4: create measure gauge for masks which will count the number
            % of blobs
            % 5: 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % vertical start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1: construct the data source nozzle
            SUBimageSource = dataSource(FileList,@(X)double(imread(X)),1);
            % if not-compiled - then full load into dev machine
            if ~isdeployed()
                SUBimageSource.fullLoad();
            end
            % 1: construct the mask nozzle
            SUBmaskSource = dataSource(mFileList,@(X)imread(X),1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2: create the thumbnail image and depth stack the images
            stacked_subMasterNozzle = dataNozzle(@(X,e0,e1)func_resizeDepthStack(X,PER,false,toSort),SUBimageSource,1);
            % create data nozzle which resizes and stacks
            nestedFunc1 = @(X,e0,e1)func_resizeDepthStack(logicalNozzle(X,@(X,Y)any(X,Y),1,false),PER,false,false);
            % create vertical OR nozzle for masks
            verticalANYnozzle = dataNozzle(nestedFunc1,SUBmaskSource,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create blob counting gauge for vertical nozzle
            blobCountGauge = nozzleGauge(verticalANYnozzle,{@(X,v)countParticleBlobs(X,v)});
            % measure max blobs
            blobM = blobCountGauge.takeMetric();
            repeatingUnits = blobM{1};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obtain matching spray nozzle
            % THIS NEED TO MADE INTO AN OBJECT
            LOOP = 2;
            nGRPs = [0 0];
            nCOMP = 5;
            [tester] = shapeVerticalStripNozzle(stacked_subMasterNozzle,verticalANYnozzle,nCOMP,tPROB,repeatingUnits,nGRPs,LOOP);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % test save
            % 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % vertical end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % horizontal start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    
            %%%% HOR
            nHORU = 1;
            hPER = .25;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1: construct the data source nozzle for resize
            % this function will resize the data from the data source
            func = @(X,e0,e1)func_thumbNail(X,hPER,false,false);
            resizeNozzle = dataNozzle(func,SUBimageSource,1);
            resizeMaskNozzle = dataNozzle(func,SUBmaskSource,1);

            
            
            % create data nozzle which resizes and stacks
            nestedFunc2 = @(X,e0,e1)func_thumbNail(logicalNozzle(X,@(X,Y)any(X,Y),1,false),hPER,false,false);


            % create vertical OR nozzle for mask
            % NEED TO HAVE THIS A GAUGE FOR PERCENT OF AVERAGE WIDTH
            reSizeH = 800;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % nozzle that will apply the attached nestedFunc2
            % this function will take the ANY operator along the horizontal
            % direction, resize it and return it as the selection object 
            % for the below nozzleManifolds
            verticalANYnozzle2 = dataNozzle(nestedFunc2,SUBmaskSource,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % nozzle manifold has access to two sources of data
            % 1: resize (data/mask) nozzle and 2: selector nozzle
            % the function attached to the nozzle will find the objects
            % from the selector nozzle,crop out data from the (data/mask)
            % nozzle and then reshape the data/mask informaiton.
            vertStripImage = nozzleManifold({resizeNozzle,verticalANYnozzle2},@(X,e0,e1)cropANDshape(X,e0,e1,reSizeH),nHORU);
            vertStripMask = nozzleManifold({resizeMaskNozzle,verticalANYnozzle2},@(X,e0,e1)cropANDshape(X,e0,e1,reSizeH),nHORU);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            vertAny = dataNozzle(@(X,e0,e1)logicalNozzle(X,@(X,Y)any(X,Y),1,false),vertStripMask,1);

            
            % create blob counting gauge for vertical nozzle
            blobCountGauge = nozzleGauge(vertAny,{@(X,v)countParticleBlobs(X,v)});
            % measure max blobs
            blobM = blobCountGauge.takeMetric();
            repeatingUnits = blobM{1};
            
            
            
           
            LOOP = 3;
            nGRPs = [0 0];
            nCOMP = 5;
            [tester2] = shapeVerticalStripNozzle(vertStripImage,vertAny,nCOMP,tPROB,repeatingUnits,nGRPs,LOOP);


            
            tester.dataSource.dataSource.dataSource.MEMstore = {};
            tester2.dataSource.dataSource.sourceNozzles{2}.dataSource.MEMstore = {};
            % create temp data source
            %tmpSource = dataSource(FileList,@(X)double(imread(X)),0);
            % set the functions to have the temp data sources
            %tester.dataSource.dataSource.dataSource = tmpSource;
            %tester2.dataSource.dataSource.sourceNozzles{1}.dataSource = tmpSource;
            % create nozzle pointer
            %tester.nozzlePtr = nozzlePtr(numel(FileList),1);
            %tester3 = dataNozzle(@(X,e0,e1)func_thumbNail(X,[1 4992],false,false),tester,1);
            %tester2.dataSource.dataSource.sourceNozzles{2}.dataSource = tester3;
            %tester2.nozzlePtr = nozzlePtr(numel(FileList)*3,1);
            %tester2.dataSource.dataSource.readsPerSource = 3;
            %tester2.dataSource.dataSource.nozzlePtr = nozzlePtr(1,3);
            
            T{t}{c}.name = [cdir(t).name '__' cropboxPath(c).name];
            T{t}{c}.function{1} = tester;
            T{t}{c}.function{2} = tester2;
            clear tester tester2 SUBimageSource
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
    T = generateImageClass(rootPath,1);
%}