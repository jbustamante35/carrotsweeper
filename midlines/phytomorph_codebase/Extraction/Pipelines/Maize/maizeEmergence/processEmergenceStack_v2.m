function [] = processEmergenceStack_v2(imageStack,oPath)
    fprintf(['****************************************************************************************\n']);
    versionString = ['Starting emergence analysis algorithm. \nPublication Version 1.0 - Monday, April 3, 2017. \n'];
    fprintf(versionString);
    fprintf(['****************************************************************************************\n']);
    
    N = 10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rectify images auto on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['Starting image rectification process \n']);tm = clock;
    [rec] = getRectification(imageStack{1},1,1);
    fprintf(['Ending image rectification process:' num2str(etime(clock,tm)) '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['Starting circle finding process \n']);tm = clock;
    [circleMatrix MASK tMASK OFFSET] = getCircles(imageStack,rec,N,168);
    fprintf(['Ending circle finding process:' num2str(etime(clock,tm)) '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
     
     
    [MAP] = getMap_ver0(tMASK,rec,OFFSET);
    [R] = getBoundingBoxes(MASK);
    for e = 1:numel(R)
        delta = bsxfun(@minus,MAP.centers,R(e).Centroid);
        [J,sidx(e)] = min(sum(delta.*delta,2));
    end
    R(sidx) = R;
    isidx = 1:numel(R);
    isidx = isidx(sidx);
    %{
    I = imread(imageStack{2});
    imshow(I,[])
    %{
    hold on
    for e = 1:numel(R)
        text(R(e).Centroid(1),R(e).Centroid(2),num2str(e))
    end
    %}
    
    hold on
    plot(MAP.centers(:,1),MAP.centers(:,2),'r*')
    for l = 1:numel(MAP.labels)
        text(MAP.centers(l,1),MAP.centers(l,2),[num2str(l) '-' MAP.labels{l}])
    end
    
    % load quick
    binarySig = csvread('/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/20170220_Camera4.csv');
    %}
    
    para.scales.value = 9;
    para.resize.value = .75;
    BORDER = 100;
    
    clear binarySig rawSignal thresh
    binarySig = zeros(numel(imageStack)-1,numel(R));
    rawSignal = binarySig;
    thresh = zeros(size(numel(R)));
    [pth,nm,ext] = fileparts(imageStack{1});
    fidx = strfind(pth,filesep);
    parfor e = 1:numel(R)
        tm = clock;
        tmpBB = R(e).BoundingBox;
        [miniStack miniMask] = diskCrop(imageStack,MASK,tmpBB,BORDER,rec);
        % measureMovie(miniStack,miniMask,0)
        % viewData(outFile,[],1,0,0);
        outFile = [oPath filesep pth((fidx(end)+1):end) '_' num2str(e) '.mat'];
        spoolToDisk_emergence(outFile,miniStack,miniMask);
        
        %squareMini(miniStack,miniMask)
        %save([oPath num2str(e)
        %[binarySig(:,e),rawSignal(:,e),thresh(e)] = processMiniStack(miniStack,miniMask,30,5,50);
        %[percentRed(:,e),ampSig(:,e)] = processMiniStack(miniStack,miniMask,30,5,50);
        %{
        %[featureStack] = surKurStack(miniStack,miniMask,para);
        %[kidx(e) Z{e}] = extractPointsFromFeatureStack(featureStack,miniMask,miniStack);
        %showAndPlot(miniStack,squeeze(miniStack(:,:,1,:)),miniMask,sK(:,1,e))
        %[kidx(e) L(e,:) bc(e,:) di(e,:)] = findGermFrame_ver0(miniStack,miniMask);
        %}
        fprintf(['Time total estimate:' num2str(numel(R)*etime(clock,tm)/12/60) '\n'])
    end
    
    
    
    %% generate graphic for presentation
    %{
    sel = 95;
    displayMaskNumbers(MASK,imageStack{3});
    tmp = rawSignal(1:end,sel);
    figure;
    plot(tmp,'k')
    hold on
    plot(thresh(sel)*ones(size(tmp)),'r')
    figure;
    fidx = find(binarySig(:,sel));
    if ~isempty(fidx)
        frame2 = fidx(1);
        F1 = miniStack(:,:,:,1);
        FE = miniStack(:,:,:,end);
        FPTPOP = miniStack(:,:,:,fidx(1)+1);
        FPOP = miniStack(:,:,:,fidx(1));
        FPPOP = miniStack(:,:,:,fidx(1)+2);

        TOT = cat(2,F1,FPOP,FPTPOP,FPPOP,FE);
        imshow(TOT/255,[]);
    end
    
    fidx = find(binarySig(:,e));
    
    for e = 1:size(miniStack,4)
        imshow(miniStack(:,:,:,e)/255,[]);
        title(num2str(e))
        drawnow
        pause(.1)
        %if e == fidx(1)
        %    waitforbuttonpress
        %end
    end
    
    g = cat(2,MOV(:,:,:,1),MOV(:,:,:,45),MOV(:,:,:,50),MOV(:,:,:,130));
    imshow(g,[]);
    imwrite(g,'/mnt/spaldingdata/nate/f1.tif');
    %}
    %{
    [pth,nm,ext] = fileparts(imageStack{1});
    fidx = strfind(pth,filesep);
    
    outFile = [oPath filesep pth((fidx(end)+1):end) '_percentRed.csv'];
    csvwrite(outFile,percentRed);
    
    outFile = [oPath filesep pth((fidx(end)+1):end) '_stdColor.csv'];
    csvwrite(outFile,ampSig);
    %}
    
    
    %% generate graphs
    genG = 0;
    if genG
        ipdf = gradient(imfilter(mean(binarySig,2),fspecial('average',[5 1]),'replicate'));   
        ipdf = ipdf / 30  * 60 ;
        ipdf(ipdf < 0) = 0;
        xlab = 1:numel(ipdf);
        xlab = (xlab*30/60 + 80)/24;
        close all
        [AX,H1,H2] = plotyy(xlab,ipdf*100,xlab,mean(binarySig,2));
        
        set(get(AX(1),'Ylabel'),'String','percent germ per hour') 
        set(get(AX(2),'Ylabel'),'String','total percent germ')
        xlabel('days after planting');
    end
    
    
    
    %generate movie
    genM = 0;
    if genM
        MOV = generateMovie(imageStack,MASK,binarySig,sidx);
        %MOV = miniStack/255;
        close all
        figure;
        nipdf = bindVec(ipdf)*size(MOV,1);
        nX = linspace(1,size(MOV,2),numel(ipdf));
        nW = bindVec(mean(binarySig,2))*size(MOV,1);
        writerObj = VideoWriter('/mnt/spaldingdata/nate/forTalk_ver1.avi');
        writerObj.FrameRate = 15;
        axis tight
        open(writerObj);
        for loop = 1:1
            for e = 1:size(MOV,4)
                imshow(MOV(:,:,:,e),[])
                title(num2str(e))
                %pause(1)
                frame = getframe;
                writeVideo(writerObj,frame);
                %{
                hold on


                %[AX,H1,H2] = plotyy(nX,nipdf,nX,nW);
                AX1 = axes;
                H1 = plot(AX1,nX,nipdf,'LineWidth',4);
                AX2 = axes
                H2 = plot(AX2,nX,nW,'LineWidth',4);            
                set(AX1,'Color','None')
                set(AX2,'Color','None')

                %}
                 
                drawnow
            end
            close(writerObj);
        end
    end
end
    %{
    

   
    
    
    
    
    
    
    %{
    toSNIP = 2;
    [tmpL] = snipSignal(L,toSNIP,7);
    [tmpbc] = snipSignal(bc,toSNIP,7);
    [tmpdi] = snipSignal(di,toSNIP,7);
    
    tmpS = tmpL.*tmpdi;
    [pidx pv] = findSignalPeak(tmpS,20);
    [bk] = getSigalBackground(tmpS,30);
    
    nsig = (pv.*bk.^-1);
    kidx = kmeans(nsig,2);
    
    for e = 1:numel(kidx)
        framesToGet = numel(imageStack);
        [frameCrop] = diskCropAtFrame(imageStack,framesToGet,R(e).BoundingBox);
        imshow(frameCrop,[]);
        title(num2str(kidx(e)))
        waitforbuttonpress
    end
    
    
    
    [bS bC bU bE bL bERR bLAM] = PCA_FIT_FULL(tmpbc,3);
    toProj = tmpS;
    oL = bsxfun(@minus,toProj,mean(toProj,1));
    oL = bsxfun(@plus,oL,bU);
    [hC] = PCA_REPROJ(oL,bE,bU);
    mL = PCA_BKPROJ(hC,bE,bU);
    mL = bsxfun(@minus,mL,bU);
    mL = bsxfun(@plus,mL,mean(toProj,1));
    
    nL = -(toProj - mL);
    
    %{
    for e = 1:numel(R)
        [bc(e,:)] = getBrightnessCurve(mov)
    end
    %}
    
    
    % feed orginal signal
    nL = L;
    tmpbc = bc;
    
    
    % feed tmp signal
    nL = tmpL;
    tmpbc = tmpbc;
    
    % mean by std
    nL = tmpL.*tmpdi;
    
    % percent threshold
    perTHRESH = .3;
    
    %nL = tmpS;
    
    % from curvature
    nL = squeeze(sK(:,1,:))';
    
    kidx = [];
    sL = [];
    for e = 1:size(nL,1)
        %[kidx(e) sL(e,:)] = find_L_threshold(nL(e,:),tmpbc(e,:),perTHRESH,1);
        [kidx(e) sL(e,:)] = find_L_threshold(nL(e,:),[],perTHRESH,1);
    end
    %}
    
    %kidx = pidx;
    %sL = tmpS;
    %OFFSET = 16;
    kidx(kidx~=0) = kidx(kidx~=0) + 0 + OFFSET;
    kidx = round(kidx);
    for e = 1:size(nL,1)
        tmpStack = [];
        v = max(sL(e,:));
        z = zeros(size(sL(e,:)));
        if kidx(e) ~= 0
            z(kidx(e)) = v;
            framesToGet = [(kidx(e)-1) kidx(e) (kidx(e)+1) numel(imageStack)];
            idx = (framesToGet > numel(imageStack));
            framesToGet(idx) = numel(imageStack);
            [frameCrop] = diskCropAtFrame(imageStack,framesToGet,R(e).BoundingBox,BORDER);
            frameCrop(:,:,:,end) = 255*flattenMaskOverlay(double(frameCrop(:,:,:,end))/255, logical(Z{e}),.25,'r');
            for i = 1:size(frameCrop,4)
                tmpStack = cat(2,tmpStack,frameCrop(:,:,:,i));
            end
        else
            framesToGet = [numel(imageStack)];
            [tmpStack] = diskCropAtFrame(imageStack,framesToGet,R(e).BoundingBox,BORDER);
            tmpStack = flattenMaskOverlay(double(tmpStack)/255, logical(Z{e}),.25,'r');
        end
        figure
        plot(sL(e,:),'k');
        hold on
        plot(z,'r')
        figure
        imshow(tmpStack);
        title(num2str(framesToGet));
        waitforbuttonpress
        hold off
        close all
        
    end
    
    
    
    
end

%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) scan for images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FilePath = '/home/nate/Downloads/Overhead_Compilation/';
    FilePath = '/home/nate/Downloads/Angle_Compilation/';
    FilePath = '/home/nate/Downloads/emergance/';
    FilePath = '/home/nate/Downloads/20151222_Camera1/';
    FilePath = '/mnt/scratch1/phytoM/flashProjects/workWithGustin/20160106_Camera1/';
    FilePath = '/mnt/scratch1/phytoM/flashProjects/workWithGustin/Checkerboard/';
    FilePath = '/mnt/scratch1/phytoM/flashProjects/workWithGustin/20160115_Camera1/';
    FilePath = '/mnt/snapper/nate/forEmergance/20160928_Camera3/';
    FilePath = '/mnt/snapper/nate/forEmergance/20161110_Camera3/';
    FilePath = '/mnt/snapper/nate/forEmergance/20161202_Camera3/';
    FilePath = '/home/nate/Downloads/20161219_Camera3/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170131_Camera3/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170131_Camera4/';

    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera1/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera2/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera3/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera4/';


    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera1/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera2/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera3/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera4/';

    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera1/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera2/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera3/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera4/';
    FileList = {};
    FileExt = {'tiff','TIF','tif','JPG','jpg'};
    FileExt = {'tiff'};
    verbose = 1;
    FileList = gdig(FilePath,FileList,FileExt,verbose);
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 0) viewDrift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = {};
    for e = 1:numel(FileList)
        [p n{e} ex] = fileparts(FileList{e});
        n{e} = str2num(n{e});
    end
    [n sidx] = sort(cell2mat(n));
    FileList = FileList(sidx);

    tmp = imread(FileList{2});
    [J,BOX] = imcrop(tmp);
    ROW = [BOX(2) BOX(2) + BOX(4)];
    COL = [BOX(1) BOX(1) + BOX(3)];

    IJ = [];
    for e = 1:numel(FileList)
        tmp = imread(FileList{e},'PixelRegion',{round(ROW) round(COL)});
        %tmp = imcrop(tmp,BOX);
        IJ(:,:,:,e) =tmp;
        imshow(tmp,[]);
        title(num2str(e));
        drawnow
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) checkerboad check
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %checkerBoardChecker(FileList);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) sort the images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oPath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
    n = {};
    for e = 1:numel(FileList)
        [p n{e} ex] = fileparts(FileList{e});
        n{e} = str2num(n{e});
    end
    [n sidx] = sort(cell2mat(n));
    FileList = FileList(sidx);
    FileList(204:end) = [];
    processEmergenceStack_v2(FileList,oPath);










%}
    %}