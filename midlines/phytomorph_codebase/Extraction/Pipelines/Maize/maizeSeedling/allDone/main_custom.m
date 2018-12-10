%% local focus
FilePath = '/mnt/tetra/nate/seedlingDATApile/maizeData/seedlingData/';
FileList = {};
FileExt = {'nef'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auto clicks for RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pad = 50;
QRSHEET = [];
QR_PIECES = [];

cnt = 1;
clickedFileList = {};
clickPoints = {};

for e = 1:1:numel(FileList)
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
        %clickedFileList{cnt} = imread(FileList{e});
        QR{cnt} = qr;
        cnt = cnt + 1;

        hold off;
        drawnow
    catch
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make noise images for red sparkle - ver2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
SPboxSequence = [1700 1700];
SPzoomSequence = [.15];

cnt = 1;
NEWY = [];
NEWX = [];
MAHGY = [];
CENTER = [];
for e = 1:100%numel(QR)
    for m = 1:1
        mag = (.8*rand(1)+.6);
        mag = 1;
        
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = double(QR{e});
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        
        
        nIMGx = rand([fliplr(SPboxSequence) 3]);
        nIMGx = imresize(nIMGx,SPzoomSequence);
        nIMGx = rand(size(nIMGx));
        
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k));
        end
        
        
        GOOD_red_corners = clickPoints{e}*SPzoomSequence_tmp;
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        
        
        
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        drawnow
        hold off
        
        
        
        
        
        
        NEWX(:,:,:,cnt) = nIMGx;
        MAHGY(cnt) = mag;
        CENTER(:,cnt) = mean(GOOD_red_corners,1);
        NEWY(:,:,cnt) = GOOD_red_corners;
        
        
        
        
        
        
        %GOOD_red_corners = bsxfun(@minus,GOOD_red_corners,mean(GOOD_red_corners,1));
        %NEWY(cnt,:) = GOOD_red_corners(:);%*SPzoomSequence^-1;
        cnt = cnt + 1;
    end
    e
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make noise images for red sparkle - ver2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
SPboxSequence = [1300 1300];
SPzoomSequence = [.15];
cnt = 1;
NEWY = [];
NEWX = [];
window_size = [15 15 3];
dilateAmount = 3;
sampleN_0 = 200;



for e = 1:100%numel(QR)
    
    for m = 1:5
        mag = (.8*rand(1)+.6);
        mag = 1;
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = double(QR{e});
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        
        
        nIMGx = rand([fliplr(SPboxSequence) 3]);
        nIMGx = imresize(nIMGx,SPzoomSequence);
        
        
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k));
        end
        
        
        GOOD_red_corners = clickPoints{e}*SPzoomSequence_tmp;
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        
        
        
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        drawnow
        hold off
        
        
        tmpF = im2colF(nIMGx,window_size,[1 1 1]);
        
        class_label = [];
        
        for class = 1:size(GOOD_red_corners,1)
            tmp_mask = zeros(size(nIMGx,1),size(nIMGx,2));
            tmp_mask(round(GOOD_red_corners(class,2)),round(GOOD_red_corners(class,1))) = 1;
            tmp_mask = imdilate(tmp_mask,strel('disk',dilateAmount,0));
            tmp_mask = im2colF(tmp_mask,window_size(1:2),[1 1]);
            tmp_mask = tmp_mask((end-1)/2,:);
            class_label(class,:) = tmp_mask;
            
            
            fidx1 = find(tmp_mask==1);
            fidx0 = find(tmp_mask==0);
            fidx0 = fidx0(randperm(numel(fidx0)));
            
            
            samp1 = tmpF(:,fidx1);
            samp0 = tmpF(:,fidx0);
            
            
            if e == 1
                Xdata{class} = [];
                Ydata{class} = [];
            end
            
            Xdata{class} = [Xdata{class},samp1,samp0(:,1:sampleN_0)];
            Ydata{class} = [Ydata{class},ones(1,numel(fidx1)),zeros(1,sampleN_0)];
            
            
        end
        
        
        %NEWX(:,:,:,cnt) = nIMGx;
        %NEWY(:,:,cnt) = GOOD_red_corners;
        %GOOD_red_corners = bsxfun(@minus,GOOD_red_corners,mean(GOOD_red_corners,1));
        %NEWY(cnt,:) = GOOD_red_corners(:);%*SPzoomSequence^-1;
        cnt = cnt + 1;
    end
    e
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for class = 1:numel(Xdata)
    sz = size(Xdata{class});
    Xdata{class} = reshape(Xdata{class},[window_size sz(end)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% class approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
net = {};
for class = 1:numel(Xdata)
    close all
    options = trainingOptions('sgdm',...
        'InitialLearnRate',.01,...
        'MaxEpochs',10,...
        'Plots','training-progress',...
        'ExecutionEnvironment','parallel',...
        'Momentum',.9);
    layers = [
        imageInputLayer(window_size,'Normalization','None')
        convolution2dLayer(5,8,'Padding','same')
        reluLayer()
        fullyConnectedLayer(2)
        softmaxLayer()
        classificationLayer()
        ];
    net{class} = trainNetwork(Xdata{class},categorical(Ydata{class}'),layers,options);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% class approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
for e = 1000%numel(QR)
    
        mag = (.8*rand(1)+.6);
        mag = 1;
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = double(QR{e});
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        
        nIMGx = rand([fliplr(SPboxSequence) 3]);
        nIMGx = imresize(nIMGx,SPzoomSequence);
        
        
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k));
        end
        
        GOOD_red_corners = clickPoints{e}*SPzoomSequence_tmp;
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        
        tmpF = im2colF(nIMGx,window_size,[1 1 1]);
        tmpF = reshape(tmpF,[window_size size(tmpF,2)]);
        M = [];
        for class = 1:numel(net)
            preY = net{class}.predict(tmpF);
            preY = preY(:,2);
            [~,midx] = max(preY);
            preY = col2im(preY,window_size(1:2),[size(nIMGx,1) size(nIMGx,2)]);
            [M(class,2),M(class,1)] = ind2sub(size(preY),midx);
           
        end
        M = bsxfun(@plus,M,(window_size(1:2)-1)/2);
        
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        plot(M(:,1),M(:,2),'go')
        drawnow
        hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train regress center of mass and resize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
options = trainingOptions('sgdm',...
    'InitialLearnRate',.0000001,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel',...
    'Momentum',.9);
szX = size(NEWX);
layers = [
    imageInputLayer(szX(1:2),'Normalization','None')
    convolution2dLayer(21,5)
    reluLayer()
    fullyConnectedLayer(3)
    regressionLayer()];
Y = [CENTER' MAHGY'];
%[nY,uY,sY] =  zscore(Y);
MAG_LOC_net = trainNetwork(NEWX(:,:,1,:),Y,layers,options);
%MAG_LOC_net = trainNetwork(NEWX(:,:,1,:),Y,MAG_LOC_net.Layers,options);
%%
%nnY = bsxfun(@minus,Y,uY);
%nnY = bsxfun(@times,nnY,sY.^-1);
MAG_LOC_net.predict(NEWX(:,:,1,1))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate grid for scatter zoom for pots
rawI = double(imread(FileList{e}))/255;
%%
rI = imresize(rawI,.15);
%rI = imresize(rI,.7);
%%
close all

[potY,potSZ] = genIXgrid2(size(rI),[size(NEWX,1)/4 size(NEWX,2)/4],[0 0],[(size(NEWX,1)-1)/2 (size(NEWX,2)-1)/2]);

imshow(rI,[]);hold on;
plot(potY(:,1),potY(:,2),'g*')
close all
N = 5;
boxSequence = repmat(boxSequence,[N 1]);
%potY = flipdim(potY,2);
[potY] = runZoomSequence2(rI,MAG_LOC_net,potY,boxSequence,ones(1,N),uY,sY,ones(1,N));
potY = flipdim(potY,2);
close all
imshow(rI,[]);
hold on
for path = 1:size(potY,1)
    plot(squeeze(potY(path,1,:)),squeeze(potY(path,2,:)),'k')
end
plot(potY(:,1,1),potY(:,2,1),'g*')
plot(potY(:,1,end),potY(:,2,end),'r*')


