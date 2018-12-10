%% load video
videoFile = '/mnt/spaldingdata/nate/mirror_images/maizeData/flights/2/Arlington/Flight 5.MP4';
videoFile = '/home/nate/Downloads/vid.m4v';
videoFile = '/home/nate/Downloads/DJI_0002.MP4';
videoFile = '/mnt/snapper/nate/droneK/20170515_17ASH_LF_AM_15m_9kmh_x5r_video_-70/DJI_A01733_C001_20170515.MP4';
%video = mmread(videoFile,1);
readerObj = VideoReader(videoFile);
%% scan for images
FilePath = '/mnt/snapper/nate/droneK/20170515_17ASH_LF_AM_15m_9kmh_x5r_video_-70/01/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
for e = 1:numel(FileList)
    [p,n,ext] = fileparts(FileList{e});
    nm(e) = str2num(n);
end
[nm,sidx] = sort(nm);
FileList=  FileList(sidx);
%%
r = {};
for e = 1:numel(FileList)
    I = imread(FileList{e});
    [r{e} c{e} v{e}] = impixel(I);
    %imshow(I,[]);
    drawnow
end
%%
for e = 1:numel(v)
    if ~isempty(v{e})
        kp(e) = 1;
    else
        kp(e) = 0;
    end
end
%%
FileList(kp==0) = [];
v(kp==0) = [];
FileList = FileList(1:numel(v));
%% make theta
close all
I = imread(FileList{1});
R = min(size(I(:,:,1)));
R = ((R^2)/2)^.5/2;
cnt = 1;
STACK = [];
for e = 5:(numel(v)-5)
    tmp = [diff(r{e},1,1) diff(c{e},1,1)];
    th(e) = atan2(-tmp(2),tmp(1))*180/pi;
    I = imread(FileList{e});
    sz = fliplr(size(I(:,:,1)));
    BOX = [sz/2 - R 2*R 2*R];
    I = imrotate(I,-th(e),'crop');
    RO = linspace(-pi,pi,20);
    for ro = 1:numel(RO)
        t = imrotate(I,RO(ro));
        STACK(:,:,:,cnt) = imresize(imcrop(t,BOX),.05);
        ANG(cnt) = RO(ro);
        cnt = cnt + 1;
        %{
        imshow(t,[]);
        hold on
        rectangle('Position',BOX,'EdgeColor','r');
        hold off
        drawnow
        %}
    end
    e
end
%%
layers = [ ...
    imageInputLayer([size(STACK,1) size(STACK,2) 3])
    convolution2dLayer(20,25)
    reluLayer
    convolution2dLayer(5,25)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','MaxEpochs',3,'InitialLearnRate', 1e-6);
net = trainNetwork(STACK/255,ANG',layers,options);
%%
pp = predict(net,STACK/255);

%% scan for MP4
FilePath = '/mnt/snapper/nate/droneK/20170515_17ASH_LF_AM_15m_9kmh_x5r_video_-70//';
FileList = {};
FileExt = {'MP4'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% read in many videos
parfor e = 1:numel(FileList)
   D{e} = readVideoMP4(FileList{e},.25);
end
%%
miniStack = [];
for e = 1:numel(D)
    miniStack = cat(4,miniStack,D{e});
end
%% make GMMM
SKIP = 20;
[GMModel] = makeGMMset(miniStack(:,:,:,1:SKIP:end));
%% scan for MP4
FilePath = '/mnt/snapper/nate/droneK/20170515_17ASH_LF_AM_15m_9kmh_x5r_video_-70/01/';
FileListT = {};
FileExt = {'tif'};
FileListT = gdig(FilePath,FileListT,FileExt,1);

%%
cI = imread(FileListT{1});
cI = zeros([size(cI) numel(FileListT)]);
NL = zeros([size(cI,1) size(cI,2) 4 numel(FileListT)]);
parfor e = 1:numel(FileListT)
    [cI(:,:,:,e),NL(:,:,:,e)] = clusterImage(FileListT{e},GMM);
    e
end
%% cluster movie and show
[CM NL] = clusterMovie(FileList{1},GMModel,1,0);
%% overlaycluster
readerObj = VideoReader(videoFile);
cnt = 1;
SKIP = 10;
close all
for e = 100:SKIP:numel(CM)
    overlay = CM{e};
    tmp = imresize(D{1}(:,:,:,e),4);
    tmp = tmp - min(tmp(:));
    tmp = tmp / max(tmp(:));
    msk = logical(full(overlay));
    msk = logical(msk | msk ~= ~bwareaopen(~msk,5000));
   % D{1}i
   
    out = flattenMaskOverlay(tmp,msk);
    imshow(out,[]);
    drawnow
    cnt = cnt + 1;
end
%% read from video object
for e = 1:1000
%while readerObj.hasFrame
    video.frames(e).cdata = readerObj.readFrame;
    e
end
%% view from video read
close all
SKIP = 10;
for e = 1:SKIP:numel(video.frames)
    imshow(double(video.frames(e).cdata)/255,[]);
    drawnow
end
%% stack video

%% get file list 
FilePath = '/mnt/snapper/nate/mirror_images/wheat/vid2/';
FilePath = '/mnt/snapper/nate/mirror_images/wheat/vid3/';
FileList = {};
FileExt = {'tif'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%  sort
for e = 1:numel(FileList)
    [pth,tmp,ext] = fileparts(FileList{e});
    nm(e) = str2num(tmp);
end
[J,sidx] = sort(nm);
FileList = FileList(sidx);
%%
sz = size(stack);
G = zeros(sz(1),sz(2),1,sz(4));
for e = 1:size(stack,4)
    G(:,:,1,e) = rgb2gray(stack(:,:,:,e)/255);
    e
end
%% get crop box
[stack] = grabStack(FileList,100,1,101);
[J BOX] = imcrop(stack(:,:,:,1)/255);
%% get stack
[stack] = grabStack(FileList,2000,1,2000+20,BOX);
%% try sheets sample
close all
%[dX dY V] = impixel(stack(:,:,:,1)/255);
RAD = 350;
[n2 n1] = ndgrid(linspace(0,RAD,100),linspace(0,pi,100));
X = n2.*cos(n1) + size(stack,2)/2;
Y = n2.*sin(n1) + size(stack,1)/2;
sheet = cat(3,X,Y);
imshow(stack(:,:,:,1)/255);
hold on
for e = 1:size(X,2)
    plot(sheet(:,e,1),sheet(:,e,2),'r')
end
%% view movie
close all
for e = 1:size(G,4)
    imshow(imresize(G(:,:,:,e),.25),[])
    drawnow
end
%%
figure
imshow(squeeze(G(:,:,:,1:95:end)),[])
%%
%%waitforbuttonpress
close all
[m g f] = lineSample(G,sheet);
imshow(stack(:,:,:,1)/255);
hold on
for e = 1:size(X,2)
    plot(sheet(:,e,1),sheet(:,e,2),'r')
end
plot(sheet(:,m,1),sheet(:,m,2),'c')
plot(-sheet(:,m,2),sheet(:,m,1),'c')
%%
getFlyDirection(FileList);
%% stack mini video for gmm
close all
masterStack = [];
TOT = 400;
disp = 1;
for e = 1:10:TOT
    e
    %osz = size(video.frames(e).cdata);
    I = imread(FileList{e});
    I = imfilter(I,fspecial('disk',21),'replicate');
    %I = rgb2hsv_fast(I);
    osz = size(I);
    I = imresize(I,.15);
    if disp
        imshow(I,[]);
        drawnow
    end
    %I = rgb2hsv(I);
    sz = size(I);
    %{
    G = [];
    for i = 1:size(I,3)
        [g1 g2] = gradient(I(:,:,i));
        g1 = reshape(g1,[sz(1)*sz(2) 1]);
        g2 = reshape(g2,[sz(1)*sz(2) 1]);
        g3 = (g1.^2 + g2.^2).^.5;
        G = [G g3];
    end
    %}
    
    I = double(reshape(I,[sz(1)*sz(2) sz(3)]));
    %masterStack = [masterStack;[I G]];
    masterStack = [masterStack;[I]];
    e
end
%% fit 
GRPS = 2;
OP = statset('Display','iter');
obj = gmdistribution.fit(masterStack,GRPS,'options',OP);
%% h cluster
idx = cluster(obj,masterStack);
obj2 = gmdistribution.fit(masterStack(idx==1,:),GRPS,'options',OP);
%% label to RGB
I = imread(FileList{100});
I = imfilter(Io,fspecial('disk',5),'replicate');
%I = rgb2hsv_fast(I);
sz = size(I);
I = double(reshape(I,[sz(1)*sz(2) sz(3)]));
g = [];
for i = 1:GRPS
    g(:,i) = mvnpdf(I,obj.mu(i,:),obj.Sigma(:,:,i));
    g2(:,i) = mvnpdf(I,obj2.mu(i,:),obj2.Sigma(:,:,i));
end
% need to automate the group ID
grpID = 2;
%grpID = 1;
grpID = 3;
grpID = 1;
[v, gidx] = max(g,[],2);
[v, gidx2] = max(g2,[],2);
img = reshape(gidx,sz(1:2));
img2 = reshape(gidx2,sz(1:2));
%%
close all
imshow(img==1 & img2 == 2,[])
%% apply gaussian
%myVideo = VideoWriter('testProces2.avi','Uncompressed AVI');
%myVideo.FrameRate = 60;
%open(myVideo);
oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/flights/1/outTest_count/';
%oMov = '/home/nate/Downloads/mov2.avi';
mkdir(oPath)
close all;
SCALE = 1;
RowNumber = 6;
startP = 2000;
flag = 0;
RowNumber = 0;

BOX = {};
TOTF = 570;
TOTF = 1000;
TOTF = 200;
TOTM = clock;
disp = 0;
close all
MASK = {};
parfor e = 1:TOTF
    tm = clock;
    %tmp = mmread(videoFile,e);
    %Io = tmp.frames(1).cdata;
    %osz = size(video.frames(e).cdata);
    %I = imresize(video.frames(e).cdata,.25);
    %Io = video.frames(e).cdata;
    %Io = imresize(Io,SCALE);
    %Io = read(readerObj,e);
    Io = imread(FileList{e});
    I = imfilter(Io,fspecial('disk',21),'replicate');
    %I = rgb2hsv_fast(I);
    sz = size(I);
    I = double(reshape(I,[sz(1)*sz(2) sz(3)]));
    g = [];
    g2 = [];
    for i = 1:GRPS
        g(:,i) = mvnpdf(I,obj.mu(i,:),obj.Sigma(:,:,i));
        g2(:,i) = mvnpdf(I,obj2.mu(i,:),obj2.Sigma(:,:,i));
    end
    
    
    % need to automate the group ID
    grpID = 2;
    %grpID = 1;
    grpID = 3;
    grpID = 1;
    %g = g(:,2) + g(:,3);
    g = g(:,1).*g2(:,2);
    
    [v, gidx] = max(g,[],2);
    %[v, gidx2] = max(g2,[],2);
    %img = reshape(gidx,sz(1:2));
    %img2 = reshape(gidx2,sz(1:2));
    
    %img = reshape(v,sz(1:2));
    pdf = reshape(g,sz(1:2));
    
    
    
    [TAN,NOR,ANG,MASK,BOX{e},CP{e},isPLOT(e)] = plotID(pdf,Io,0);
    
    
    
    
    if disp
        imshow(Io,[]);
        if ~isempty(BOX{e})
            hold on
            plot(BOX{e}(:,2),BOX{e}(:,1),'r','LineWidth',3)
            hold off
            
        end
        drawnow
    end
    
    fprintf(['Done with:' num2str(e) ':' num2str(TOTF) ':' num2str(etime(clock,tm)) '\n'])
end
etime(clock,TOTM)

%%
%%
writerObj = VideoWriter('/mnt/spaldingdata/nate/mov3.avi');
writerObj.FrameRate = 10;
%axis tight
open(writerObj);
plotNumber = 1;
flagStart = 0;close all;
for e = 1:TOTF
    %tmp = mmread(videoFile,e);
    %Io = tmp.frames(1).cdata;
    Io = video.frames(e).cdata;
    %Io = flattenMaskOverlay(Io,MASK{e},.5,'r');
    imshow(Io,[]);
    
    if ~isempty(BOX{e})
        flagStart = 1;
        hold on
        plot(BOX{e}(:,2),BOX{e}(:,1),'r','LineWidth',3);
        ht = text('FontUnits','inches');
        set(ht,'position',fliplr(CP{e}));
        set(ht,'String',num2str(plotNumber));
        set(ht,'FontSize',1);
        set(ht,'BackgroundColor','w');
        get(ht,'FontUnits');
        hold off
        drawnow
    else
        if flagStart
            flagStart = 0;
            plotNumber = plotNumber + 1;
        end
    end
    title(num2str(e))
    frame = getframe;
    frame.cdata = imresize(frame.cdata,[481 861]);
    writeVideo(writerObj,frame);
    %{
    gidx = gidx == grpID;
    
    
    
    
    gidx = reshape(gidx,sz(1:2));
    gidx = imopen(gidx,strel('disk',45));
    gidx = imerode(gidx,strel('disk',21));
    gidx = imfill(gidx,'holes');
    gidx = imclearborder(gidx);
    gidx = bwareaopen(gidx,10000);
    %}
    
    
    
    %{
    
    size(Io)
    if any(gidx(:))
        [RowNumber,flag] = twist2(Io,gidx,RowNumber,flag,writerObj,0);
    else
        
        imshow(Io,[])
        %axis tight
        %pause(10)
        drawnow
        
        frame = getframe;
        frame.cdata = imresize(frame.cdata,[481 861]);
        writeVideo(writerObj,frame);
        
    end
   %close all
    %}
end
 close(writerObj);
%%

    gidx = bwareaopen(gidx,round(SCALE*1000));
    
    gidx(1,:) = 0;
    gidx(end,:) = 0;
    gidx = imclearborder(gidx);
    
    
    
    
    
    
    proj = sum(gidx,2);
    midx = imextendedmax(proj,200);
    
  
    
    
    Io = imresize(Io,SCALE.^-1);
    gidx = imresize(gidx,SCALE.^-1);
    out = flattenMaskOverlay(Io, gidx,.5,'r');     
    %figure;
    h = imshow(out,'Border','tight');
    hold on;
    %plot(proj,1:size(out,1),'b');
    %plot(midx*max(proj(:)),1:size(out,1),'k')
    
    N = 100;
    M = 100;
    CONSTX = 0;
    CONSTY = 100;
    L1 = generateLine([CONSTX+1 CONSTY],[CONSTX+1 size(gidx,1)-CONSTY],N);
    %plot(L1(1,:),L1(2,:),'r','LineWidth',5);
    L2 = generateLine([size(gidx,2)-CONSTX CONSTY],[size(gidx,2)-CONSTX size(gidx,1)-CONSTY],N);
    %plot(L2(1,:),L2(2,:),'r','LineWidth',5);
    DELTA = 100;
    SKIP = 10;
    samp = double(gidx);
    S = [];
    %pl = figure;
    ci = 1;
    cj = 1;
    E = [];
    DI = -DELTA:SKIP:DELTA;
    DJ = -DELTA:SKIP:DELTA;
    for i = DI
        d1 = [0;i];
        tL1 = bsxfun(@minus,L1,d1);
        cj = 1;
        for j = DJ
            d2 = [0;j];
            tL2 = bsxfun(@minus,L2,d2);
            
            %imshow(out);
            %hold on
            for l = 1:N
                X1 = linspace(tL1(1,l),tL2(1,l),M);
                X2 = linspace(tL1(2,l),tL2(2,l),M);
                %plot(X1,X2);
                S(l,:) = ba_interp2(samp,X1,X2);
            end
            %drawnow
            sS = sum(S,2);
            E(ci,cj) = entropy(sS);
            sS = sS / sum(sS);
            part = log2(sS);
            part(isinf(part)) = 0;
            E(ci,cj) = -sum(sS.*part);
            
            %figure(pl);
            %plot(sS);
            %drawnow
            cj = cj + 1;
        end
        ci = ci + 1;
    end
    [V,idx] = min(E(:));
    [ri ci] = find(E == V);
    d1 = [0;DI(ri(1))];
    d2 = [0;DJ(ci(1))];
    tL1 = bsxfun(@minus,L1,d1);
    tL2 = bsxfun(@minus,L2,d2);
    for l = 1:N
        X1 = linspace(tL1(1,l),tL2(1,l),M);
        X2 = linspace(tL1(2,l),tL2(2,l),M);
        XB(:,:,l) = [X1;X2];
        %plot(X1,X2,'g');
        S(l,:) = ba_interp2(samp,X1,X2);
    end
    sS = sum(S,2);
    midx = imextendedmax(sS,10);
    R = regionprops(midx,'Centroid');
    MAG = 130;
    
    
    
    noRow = 1;
    for i = 1:numel(R)
        tmp = XB(:,:,round(R(i).Centroid(2)));
        V = ba_interp2(double(gidx),tmp(1,:),tmp(2,:));
        fidx = find(V>.5);
        Rr = regionprops(V>0,'PixelIdxList','Area');
        [V,midx] = max([Rr.Area]);
        
        
        
        
        
        st = tmp(:,min(fidx));
        sp = tmp(:,max(fidx));
        delta = sp - st;
        centerPoint = st + delta*.5;
        delta = delta/norm(delta);
        nor = [delta(2);-delta(1)];
        tan = [-nor(2);nor(1)];
        centerVec = centerPoint - fliplr(size(gidx))'/2;
        quiver(size(gidx,2)/2,size(gidx,1)/2,nor(1),nor(2),1000);
        polyBox = [st + MAG*nor sp + MAG*nor sp - MAG*nor st - MAG*nor st + MAG*nor];
        
        
        dis = nor'*centerVec;
        
        
        % draw count zone
        imageCenter = fliplr(size(gidx))'/2;
        minZone = 0;
        zoneWidth = 100;
        startZone = imageCenter + nor*minZone;
        endZone = imageCenter + nor*minZone + nor*zoneWidth;
        plot([startZone(1) endZone(1)],[startZone(2) endZone(2)],'r','LineWidth',5);
        
        
        
        
        BB = [[1;1],[size(gidx,2);1],[size(gidx,2);size(gidx,1)],[1;size(gidx,1)],[1;1]];
        in = inpoly(polyBox',BB');
        
        if dis > minZone & dis < (minZone + zoneWidth);
             flag = 1;
             plot(centerVec(1),centerVec(2),'k*');
             quiver(size(gidx,2)/2,size(gidx,1)/2,centerVec(1),centerVec(2),0);
             
             
        
            if all(in)
                plot(polyBox(1,:),polyBox(2,:),'g','LineWidth',5);
                plot(tmp(1,min(fidx):max(fidx)),tmp(2,min(fidx):max(fidx)),'g');
                ht = text('FontUnits','inches');
                set(ht,'position',centerPoint);
                set(ht,'String',num2str(RowNumber));
                set(ht,'FontSize',1);
                set(ht,'BackgroundColor','w');
                get(ht,'FontUnits')
            end
            noRow = 0;
        else
           
        end
        
    end
    
    if noRow
        if flag
            RowNumber = RowNumber + 1;
            flag = 0;
        end
    end
   
    
    
    
    %print(gcf,[oPath num2str(e-startP) '.tif'],'-dtiff','-r600');
    saveas(gcf,[oPath num2str(e-startP) '.tif']);
    %frame = getframe;
    %out = imresize(out,1);
    %imwrite(out,[oPath num2str(e-1) '.png']);
    %frame = im2frame(out);
    %writeVideo(myVideo,frame);
    
    %{
    figure;
    gidx = label2rgb(gidx);
    imshow(gidx,[]);
    %}
    %figure;
    %imshow(video.frames(e).cdata,[])
    
    
    %{
    %gidx = round(imresize(gidx,osz(1:2)));
    
    drawnow
    
    
    %}
    drawnow
end
%close(myVideo);
%% apply gaussian
%myVideo = VideoWriter('testProces2.avi','Uncompressed AVI');
%myVideo.FrameRate = 60;
%open(myVideo);
oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/flights/1/outTest_count/';

mkdir(oPath)
close all;
SCALE = 1;
RowNumber = 6;
startP = 2000;
for e = 1:570%startP:1:19000%490:numel(video.frames)
    e
    osz = size(video.frames(e).cdata);
    %I = imresize(video.frames(e).cdata,.25);
    Io = video.frames(e).cdata;
    Io = imresize(Io,SCALE);
    I = imfilter(Io,fspecial('disk',3));
    %I = rgb2hsv(I);
    sz = size(I);
    %{
    G = [];
    for i = 1:size(I,3)
        [g1 g2] = gradient(I(:,:,i));
        g1 = reshape(g1,[sz(1)*sz(2) 1]);
        g2 = reshape(g2,[sz(1)*sz(2) 1]);
        g3 = (g1.^2 + g2.^2).^.5;
        G = [G g3];
    end
    %}
    
    I = double(reshape(I,[sz(1)*sz(2) sz(3)]));
    %I = [I G];
    
    
    g = [];
    for i = 1:GRPS
        g(:,i) = mvnpdf(I,obj.mu(i,:),obj.Sigma(:,:,i));
    end
    
    
    %{
    for i = 1:GRPS
        G = reshape(g(:,i),sz(1:2));
        imshow(G,[]);
        drawnow
        waitforbuttonpress
    end
    %}
    
 
    
    [v, gidx] = max(g,[],2);
    %gidx = gidx== 2 | gidx == 3;
    %gidx = gidx== 1 | gidx == 3 | gidx == 4 & v > 10^-6;
    gidx = gidx == 3 | gidx == 6 | gidx == 5;
    gidx = reshape(gidx,sz(1:2));    
    gidx = bwareaopen(gidx,round(SCALE*1000));
    
    gidx(1,:) = 0;
    gidx(end,:) = 0;
    gidx = imclearborder(gidx);
    
    
    
    
    
    
    proj = sum(gidx,2);
    midx = imextendedmax(proj,200);
    
  
    
    
    Io = imresize(Io,SCALE.^-1);
    gidx = imresize(gidx,SCALE.^-1);
    out = flattenMaskOverlay(Io, gidx,.5,'r');     
    %figure;
    h = imshow(out,'Border','tight');
    hold on;
    %plot(proj,1:size(out,1),'b');
    %plot(midx*max(proj(:)),1:size(out,1),'k')
    
    N = 100;
    M = 100;
    CONSTX = 0;
    CONSTY = 100;
    L1 = generateLine([CONSTX+1 CONSTY],[CONSTX+1 size(gidx,1)-CONSTY],N);
    %plot(L1(1,:),L1(2,:),'r','LineWidth',5);
    L2 = generateLine([size(gidx,2)-CONSTX CONSTY],[size(gidx,2)-CONSTX size(gidx,1)-CONSTY],N);
    %plot(L2(1,:),L2(2,:),'r','LineWidth',5);
    DELTA = 100;
    SKIP = 10;
    samp = double(gidx);
    S = [];
    %pl = figure;
    ci = 1;
    cj = 1;
    E = [];
    DI = -DELTA:SKIP:DELTA;
    DJ = -DELTA:SKIP:DELTA;
    for i = DI
        d1 = [0;i];
        tL1 = bsxfun(@minus,L1,d1);
        cj = 1;
        for j = DJ
            d2 = [0;j];
            tL2 = bsxfun(@minus,L2,d2);
            
            %imshow(out);
            %hold on
            for l = 1:N
                X1 = linspace(tL1(1,l),tL2(1,l),M);
                X2 = linspace(tL1(2,l),tL2(2,l),M);
                %plot(X1,X2);
                S(l,:) = ba_interp2(samp,X1,X2);
            end
            %drawnow
            sS = sum(S,2);
            E(ci,cj) = entropy(sS);
            sS = sS / sum(sS);
            part = log2(sS);
            part(isinf(part)) = 0;
            E(ci,cj) = -sum(sS.*part);
            
            %figure(pl);
            %plot(sS);
            %drawnow
            cj = cj + 1;
        end
        ci = ci + 1;
    end
    [V,idx] = min(E(:));
    [ri ci] = find(E == V);
    d1 = [0;DI(ri(1))];
    d2 = [0;DJ(ci(1))];
    tL1 = bsxfun(@minus,L1,d1);
    tL2 = bsxfun(@minus,L2,d2);
    for l = 1:N
        X1 = linspace(tL1(1,l),tL2(1,l),M);
        X2 = linspace(tL1(2,l),tL2(2,l),M);
        XB(:,:,l) = [X1;X2];
        %plot(X1,X2,'g');
        S(l,:) = ba_interp2(samp,X1,X2);
    end
    sS = sum(S,2);
    midx = imextendedmax(sS,10);
    R = regionprops(midx,'Centroid');
    MAG = 130;
    
    
    
    noRow = 1;
    for i = 1:numel(R)
        tmp = XB(:,:,round(R(i).Centroid(2)));
        V = ba_interp2(double(gidx),tmp(1,:),tmp(2,:));
        fidx = find(V>.5);
        Rr = regionprops(V>0,'PixelIdxList','Area');
        [V,midx] = max([Rr.Area]);
        
        
        
        
        
        st = tmp(:,min(fidx));
        sp = tmp(:,max(fidx));
        delta = sp - st;
        centerPoint = st + delta*.5;
        delta = delta/norm(delta);
        nor = [delta(2);-delta(1)];
        tan = [-nor(2);nor(1)];
        centerVec = centerPoint - fliplr(size(gidx))'/2;
        quiver(size(gidx,2)/2,size(gidx,1)/2,nor(1),nor(2),1000);
        polyBox = [st + MAG*nor sp + MAG*nor sp - MAG*nor st - MAG*nor st + MAG*nor];
        
        
        dis = nor'*centerVec;
        
        
        % draw count zone
        imageCenter = fliplr(size(gidx))'/2;
        minZone = 0;
        zoneWidth = 100;
        startZone = imageCenter + nor*minZone;
        endZone = imageCenter + nor*minZone + nor*zoneWidth;
        plot([startZone(1) endZone(1)],[startZone(2) endZone(2)],'r','LineWidth',5);
        
        
        
        
        BB = [[1;1],[size(gidx,2);1],[size(gidx,2);size(gidx,1)],[1;size(gidx,1)],[1;1]];
        in = inpoly(polyBox',BB');
        
        if dis > minZone & dis < (minZone + zoneWidth);
             flag = 1;
             plot(centerVec(1),centerVec(2),'k*');
             quiver(size(gidx,2)/2,size(gidx,1)/2,centerVec(1),centerVec(2),0);
             
             
        
            if all(in)
                plot(polyBox(1,:),polyBox(2,:),'g','LineWidth',5);
                plot(tmp(1,min(fidx):max(fidx)),tmp(2,min(fidx):max(fidx)),'g');
                ht = text('FontUnits','inches');
                set(ht,'position',centerPoint);
                set(ht,'String',num2str(RowNumber));
                set(ht,'FontSize',1);
                set(ht,'BackgroundColor','w');
                get(ht,'FontUnits')
            end
            noRow = 0;
        else
           
        end
        
    end
    
    if noRow
        if flag
            RowNumber = RowNumber + 1;
            flag = 0;
        end
    end
   
    
    
    
    %print(gcf,[oPath num2str(e-startP) '.tif'],'-dtiff','-r600');
    saveas(gcf,[oPath num2str(e-startP) '.tif']);
    %frame = getframe;
    %out = imresize(out,1);
    %imwrite(out,[oPath num2str(e-1) '.png']);
    %frame = im2frame(out);
    %writeVideo(myVideo,frame);
    
    %{
    figure;
    gidx = label2rgb(gidx);
    imshow(gidx,[]);
    %}
    %figure;
    %imshow(video.frames(e).cdata,[])
    
    
    %{
    %gidx = round(imresize(gidx,osz(1:2)));
    
    drawnow
    
    
    %}
    drawnow
end
%close(myVideo);