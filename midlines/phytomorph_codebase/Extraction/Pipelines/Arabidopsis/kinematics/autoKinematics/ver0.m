FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/';
%FilePath = '/home/nate/Downloads/Wassim''s compilation of old data/';
%FilePath = '/home/nate/Downloads/Wassim/';
FileList = {};
FileExt = {'tiff','TIF','png'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% 
for e = 1:numel(SET)
    SET{e}{1}
    e
end
%%
close all
[R TH] = ndgrid(linspace(min(size(I{cnt})/2),min(size(I{cnt})/2),1),linspace(-pi,pi,400));
X = R.*cos(TH) + size(I{cnt},2)/2;
Y = R.*sin(TH) + size(I{cnt},1)/2;
Z = zeros(size(I{cnt}));
Z(size(I{cnt},1)/2,size(I{cnt},2)/2) = 1;
Z = bwdist(Z);
Z = double(Z < min(size(I{cnt})/2));
alpha = pi/4*180/pi;
%W = zeros(100*100,2,numel(TH));
W = [];
for e = 1:numel(TH)
    
    %{
    [R1 TH1] = ndgrid(linspace(0,min(size(I{cnt})/2),100),linspace(-pi,pi,100));
    X1 = R1.*cos(TH1+TH(e)+pi);
    Y1 = R1.*sin(TH1+TH(e)+pi);
    %}
    
    [R1 TH1] = ndgrid(linspace(0,min(size(I{cnt})/2)*sin(alpha*pi/180),100),linspace(-(90-alpha/2)*pi/180,(90-alpha/2)*pi/180,100));
    X1 = R1.*cos(TH1+TH(e)+pi);
    Y1 = R1.*sin(TH1+TH(e)+pi);
    
    X2 = X1 + X(e);
    Y2 = Y1 + Y(e);
    
    W(:,1:2,e) = [X2(:) Y2(:)];
    %{
    V = ba_interp2(Z,X2(:),Y2(:));
    idx = find(V==1);
    W{e} = [X2(idx) Y2(idx)];
    %}
    %{
    imshow(Z,[])
    hold on
    plot(X2,Y2,'b.')
    plot(X2(end/2),Y2(end/2),'r.')
    drawnow
    hold off
    %}
end


%%
A = 50;
I = imread(SET{1}{1});
iVEC = numel(SET)-1;
iVEC = 1;%10;
iVEC = 1:numel(SET);
I = zeros([size(I) numel(iVEC)]);
I = {};
cnt = 1;
VAL = [];
close all
for e = iVEC
    I{cnt} = imread(SET{e}{1});
    
    STACK = [];
    for t = 1:20
        J = imread(SET{e}{t});
        STACK(:,:,t) = J(:,:,1);
    end
    STACK = std(double(STACK),1,3);
    WE = [];
    for w = 1:size(W,3)
        WE(:,w) = ba_interp2(STACK,W(:,2,w),W(:,1,w));
    end
    
    I{cnt} = I{cnt}(:,:,1);
    
    %{
    [R TH] = ndgrid(linspace(0,min(size(I{cnt})/2),500),linspace(-pi,pi,800));
    X = R.*cos(TH) + size(I{cnt},2)/2;
    Y = R.*sin(TH) + size(I{cnt},1)/2;
    VAL(:,:,cnt) = ba_interp2(double(I{cnt}),Y,X);
    %}
    
    
    toSAMP = imfilter(double(I{cnt}),fspecial('average',121),'replicate');
    [dx dy] = gradient(toSAMP);
    toSAMP = (dx.^2 + dy.^2).^.5;
    toSAMP = imfilter(toSAMP,fspecial('average',121),'replicate');
    if e == 1
        VALT = [];
        for w = 1:size(W,3)
            VALT(:,w) = ba_interp2(toSAMP,W(:,2,w),W(:,1,w));
        end
    end
    VAL = [];
    for w = 1:size(W,3)
        VAL(:,w) = ba_interp2(toSAMP,W(:,2,w),W(:,1,w));
    end
    
    
    DIS = zeros(size(VAL,2),1);
    for r = 1:(size(VAL,2))
        tmpV = bindVec(circshift(VAL,r-1,2));
        delta = sum(tmpV(:).*bindVec(VALT(:)));
        %delta = norm(tmpV(:)-bindVec(VALT(:)));
        DIS(r) = norm(delta(:));
        r
    end
    DIS = bindVec(DIS).*bindVec(imcomplement(mean(WE,1)))';
    [~,midx] = max(DIS);
    tmp{e} = imrotate(I{e},(midx-1)*2*pi/400*180/pi,'crop');
    
    imshow(tmp{e},[])
    drawnow
    cnt
    numel(SET)
    %imshow(I{cnt},[]);
    cnt = cnt + 1;
    %drawnow
end
%% 
toT = 1;
for e = 1:size(VAL,3)
    
    for r = 0:(size(VAL,2)-1)
        DIS(r+1,w) = 0;
        for w = 1:size(VAL,2)
            DIS(r+1,w) = norm(circshift(VAL(:,w,e),r,2) - VAL(:,w,toT));
        end
        r
    end
    
    [~,midx] = min(sum(DIS,1));
    tmp{e} = imrotate(I{e},-(midx-1)*2*pi/400,'crop');
    
end
%%
A = size(VAL,1)-1;
SKIP = 10;
SV = 21;
toT = 1;
target = VAL((end-A):SKIP:end,:,toT);
target = [target,target,target];
target = imfilter(target,fspecial('disk',SV),'replicate');
target = gradient(target);
target = target(:,(size(VAL,2)+1):(size(VAL,2)*2));
for e = 1:size(VAL,3)
    
    for r = 1:size(VAL,2)
        source = VAL((end-A):SKIP:end,:,e);
        source = circshift(source,r,2);
        source = [source,source,source];
        source = imfilter(source,fspecial('disk',SV),'replicate');
        source = gradient(source);
        source = source(:,(size(VAL,2)+1):(size(VAL,2)*2));
        l(r) = norm(target(:) - source(:));
        %imshow(source,[])
        %drawnow
    end
    e
    [~,midx] = min(l);
    VAL(:,:,e) = circshift(VAL(:,:,e),midx-1,2);
    tmp{e} = imrotate(I{e},-midx*2*pi/800,'crop');
    imshow(tmp{e},[]);
    drawnow
end
%%

%%
close all
for e = 1:size(I,3)
    imshow(I(:,:,e),[])
    drawnow
end
%%
para.resize.value = .15;
para.scales.value = [10];
[K] = surKur(I(:,:,e),para);
%% look at disk
close all
for e = 1:10

    
LOC = SZ/2;
%[LOC(2) LOC(1) V] = impixel(I(:,:,e));
SZ = [size(I,1) size(I,2)];
R = min(SZ/2);
dR = 600;
dR = R;
[n1 n2] = ndgrid(linspace(R-dR,R,600),linspace(-pi,pi,1000));
X = n1.*cos(n2) + LOC(2);
Y = n1.*sin(n2) + LOC(1);
F = ba_interp2(double(I(:,:,e)),X,Y);
close all
imshow(I(:,:,e),[]);
hold on
for s = 1:10:size(n1,2)
    plot(X(:,s),Y(:,s),'r');
end
for s = 1:10:size(n1,1)
    plot(X(s,:),Y(s,:),'r');
end
fF = [];
for f = 1:size(F,1)
    g = fspecial('gaussian',[1 20],7);
    fF(f,:) = imfilter(F(f,:),g,'circular');
    fF(f,:) = edge(F(f,:));
end
fF = imclose(fF,strel('disk',3));
fF = bwareaopen(fF,50);
fidx = find(fF);

plot(X(fidx),Y(fidx),'.')
drawnow
waitforbuttonpress
end
para.resize.value = .25;
para.scales.value = [10];
K = surKur(F,para);
figure;
imshow(fF,[])
%% get skeletons
skeleton = {};
for e = 1:numel(I)
     [skeleton{e}] = getSkeleton(I{e},0);
end
%% look at skeletons
for e = 1:numel(skeleton)
    try
        close all


        Z = zeros(size(I{e},1),size(I{e},2));
        for r = 1:4
            Z(:,1) = 1;
            Z = imrotate(Z,90);
        end
        [z1 z2] = find(Z);
        tskel = logical(skeleton{e});
        tskel = bwmorph(tskel,'thin',inf);
        tskel = bwmorph(tskel,'skeleton',inf);

        P = [];
        for k = 1:100
            ep = bwmorph(tskel,'endpoints',1);
            [e1 e2] = find(ep);
            for p = 1:2
                idx = snapTo([z1 z2],[e1(p) e2(p)]);
                dist(p) = (z1(idx) - e1(p)).^2 + (z2(idx) - e2(p)).^2;
            end
            if isempty(P)
                [P,midx] = max(dist);
                tipPoint = [e1(midx) e2(midx)];
            end
            PTS = [e1 e2];
            fidx = bsxfun(@eq,PTS,tipPoint);
            fidx = all(fidx,2);
            tipPoint = [e1(fidx) e2(fidx)];
            e1(fidx) = [];
            e2(fidx) = [];
            tskel(e1,e2) = 0;
        end


        ep = bwmorph(tskel,'endpoints',1);
        [e1 e2] = find(ep);
        PTS = [e1 e2];
        fidx = bsxfun(@eq,PTS,tipPoint);
        fidx = all(fidx,2);
        tipPoint = [e1(fidx) e2(fidx)];
        basePoint = [e1(~fidx) e2(~fidx)];


        [x,y] = find(tskel);
        % stack the skeleton points for tracing
        DP = [x y]';
        % make adjaceny matrix
        T = Radjacency(DP,3);
        % find the longest path from the stem end point to the leaf tip
        pathcost = [];
        path = {};
        [idx(1)] = snapTo(DP',tipPoint);
        [idx(2)] = snapTo(DP',basePoint);
        [path , pathcost]  = dijkstra(T , idx(2) , idx(1));
        path = DP(:,path);



        path = imfilter(path,fspecial('average',[1 100]),'replicate');
        path = arcLength(path','arcLen')';

        tskel = imdilate(tskel,strel('disk',3));
        out = flattenMaskOverlay(I{e}/255,tskel);
        imshow(out,[]);
        imshow(I{e},[]);
        hold on
        plot(path(2,:),path(1,:),'r','LineWidth',2);
        plot(path(2,1),path(1,1),'g*')
        plot(path(2,end),path(1,end),'k*')



        %{
        h1 = figure;
        h2 = figure;
        h3 = figure;
        WIND = 1:800;
        SIG = [];
        cl = [];
        cnt = 1;
        for e = 1:10:(size(FTm,1)-max(WIND))
            patch = FTm(WIND+e,:);
            [g1 g2] = gradient(patch);
            g2 = bsxfun(@minus,g2,mean(g2,1));
            SIG(:,cnt) = mean(abs(fft(g2,[],1)),2);
            figure(h1);
            imshow(patch,[]);
            figure(h2);
            imshow(SIG(1:150,:));
            minis = imfilter(SIG(:,end),fspecial('average',[5 1]),'circular');
            [J,cl(cnt)] = max(minis(1:150,end));
            figure(h3);
            plot(cl)
            drawnow
            cnt = cnt + 1;
        end
        %}

        WIDTH_NUMP = 600;
        PCA_RHO = 50;
        WIDTH = 300;
        Domain = genCurvilinearDomain(path',PCA_RHO,WIDTH,WIDTH_NUMP,[]);

        dX = -diff(Domain,1,1);
        SNIP = 200;
        dX = mean(dX(1:SNIP,:,:),1);
        dNOR = sum(dX.^2,3).^-.5;
        dX = bsxfun(@times,dX,dNOR);
        EXT = 500;
        EXT = linspace(0,EXT,EXT);
        EXT = bsxfun(@times,EXT',dX);
        EXT = bsxfun(@plus,EXT,Domain(1,:,:));


        mDomain = cat(1,flipdim(EXT,1),Domain);
        F1 = mDomain(:,:,1);
        F2 = mDomain(:,:,2);
        FT = ba_interp2(double(I{e})/255,F1,F2);
        WIDTH = mean(FT,1);
        LENGTH = mean(FT,2);
        sL = std(FT,1,2);
        thresh = graythresh(WIDTH);
        thresh2 = graythresh(sL);
        MSK = double(sL > thresh2) * double(WIDTH < thresh);
        R = regionprops(logical(MSK),'Area','PixelIdxList');
        [J,midx] = max([R.Area]);
        MSK = zeros(size(MSK));
        MSK(R(midx).PixelIdxList) = 1;


        RP = regionprops(MSK,'BoundingBox');
        FTm = imcrop(FT,RP(1).BoundingBox);
        [d1 d2] = gradient(FTm(300:end,:));
        sig = bsxfun(@minus,d1,mean((d1),2));
        k = fft((sig),[],2);
        k = mean(abs(k),1);
        pidx = find(imdilate(k,strel('disk',5)) == k);
        pidx = pidx(1)-1;
        %[v,pidx] = max(k);


        pidx = pidx - 1;
        T = size(k,2)/pidx;
        %{
        figure;
        imshow(FTm(1000:1300,:),[]);
        title(num2str(size(FTm,2)/T));
        waitforbuttonpress
        %}
        figure;
        TIP = FTm(1:500,:);
        imshow(TIP);
    
        WID = find(any(MSK,1));
       
        midline = find(any(MSK,2));
        midline = midline(1):midline(end);
        midline = sub2ind([size(mDomain,1) size(mDomain,2)],midline',round(mean(WID))*ones(numel(midline),1));
        midlineM = [];
        midlineM(:,1) = F1(midline);
        midlineM(:,2) = F2(midline);

        dB = bwboundaries(MSK);
        dB = dB{1};
        idx = sub2ind([size(mDomain,1) size(mDomain,2)],dB(:,1),dB(:,2));
        BOX = [];
        BOX(:,1) = F1(idx);
        BOX(:,2) = F2(idx);

        figure;
        imshow(FT',[]);
        hold on;
        LENGTH = 1000*(LENGTH - mean(LENGTH)) + 150;

        plot(1:numel(LENGTH),LENGTH)
        figure;
        imshow(MSK,[])
        figure;
        imshow(double(I{e})/255)
        hold on
        plotCurvilinearGrid(Domain,[30 30]);
        plotCurvilinearGrid(EXT,[30 30],[],'b');

        plot(path(2,:),path(1,:),'k','LineWidth',2);
        plot(path(2,1),path(1,1),'g*')
        plot(path(2,end),path(1,end),'b*')
        plot(BOX(:,1),BOX(:,2),'g','LineWidth',2);
        plot(midlineM(:,1),midlineM(:,2),'c','LineWidth',2);
        drawnow
        hold off

        
        waitforbuttonpress

        DST = 1;
        tSET = SET{iVEC(e)}(1:1:end);
        DSM = 30;
        diskWID = numel(WID)/2;
        para.pointList = fliplr(midlineM(1:DSM:end,:));
        para.domainPara{1}.type = 'disk';
        para.domainPara{1}.value{1} = [0 diskWID diskWID];
        para.domainPara{1}.value{2} = [-pi pi round(diskWID*2*pi)];
        para.THRESH =.0001;
        para.TIME = [];
        para.h = figure;
        [pointList] = auto_minikini(tSET,para);
        
        SAM = [1];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 2: center the tip
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Xm = pointList;        
        for t = 1:size(Xm,3)
            Xm(:,:,t) = Xm(:,:,t) - repmat(Xm(1,:,t),[size(Xm(:,:,t),1) 1 1]);
        end
        DST = 1;
        tmpX = Xm(:,:,1:DST:end);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the position values
        [J,POS,VEL] = gradient(tmpX);
        % POS = diff(tmpX,1,1); 
        % POS = gradient(tmpX);
        POS = sum(POS.*POS,2).^.5;
        h1 = mean(POS(:));
        %POS = cat(1,zeros(1,1,size(POS,3)),POS);
        POS = cumsum(POS,1);
        POS = squeeze(POS);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the velocity values
        %VEL = diff(tmpX,1,3);
        VEL = sum(VEL.*VEL,2).^.5;
        h2 = mean(VEL(:));
        VEL = squeeze(VEL);
        VEL = VEL * 1^-1;
        figure;
        plot(POS(:),VEL(:),'.')
        hold on
        TM_DOMAIN = repmat(1:size(POS,2),[size(POS,1) 1]);
        [iD1 iD2] = ndgrid(linspace(min(POS(:)),max(POS(:,1)),500),1:size(POS,2));
        
        vq = griddata(POS,TM_DOMAIN,VEL,iD1,iD2);
       
        
        %{
        POS = POS;
        VEL = VEL;
        
        
        DOMAIN = (1:size(midlineM,1));
        
        DELTA = 1;
        
        
        Zi = [];
        Zf = [];
        for tm = 1:(size(POS,2)-DELTA)
            T = tm*ones(size(DOMAIN));
            Tp = (tm+DELTA)*ones(size(DOMAIN));
            for k = 1:size(G,3)
                X = [POS(:) TM_DOMAIN(:)];
                Y = G(:,:,k);
                XI = [DOMAIN(:) T(:)];
                XI2 = [DOMAIN(:) Tp(:)];
                Zi(:,k) = griddatan(fliplr(X),Y(:),fliplr(XI),'linear');
                Zf(:,k) = griddatan(fliplr(X),Y(:),fliplr(XI2),'linear');
            end
            VEL = Zf - Zi;
            VEL = sum(VEL.*VEL,2).^.5;
            plot(DOMAIN,VEL,'k');
            drawnow
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 4: fit curves t data        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % model domain        
        tmpDOMAIN = linspace(0,3500,3500);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over each downsample value - for skipping frames
        for sam = 1:numel(SAM)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % restore whole from backup            
            tmp_Xm = Xm(:,:,1:SAM(sam):end);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop over each window value - for sliding window
            for win = 1:numel(WINDOW)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % clear vars
                totalREGR = [];
                totalVEL = [];
                totalDomain1 = [];
                totalDomain2 = [];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % loop over segments of size WINDOW(win)
                Nseg = (size(tmp_Xm,3)-WINDOW(win)-1);
                
                for seg = 1:Nseg
                    fprintf(['Done@' num2str(seg) ':' num2str(Nseg) '\n']);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % obtain segment
                    segValue = seg:seg+WINDOW(win);
                    Xwin = tmp_Xm(:,:,segValue);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % get the position values
                    POS = (Xwin(:,:,1:end-1) + Xwin(:,:,2:end))/2;
                    POS = diff(Xwin,1,1);   
                    POS = sum(POS.*POS,2).^.5;
                    POS = cat(1,zeros(1,1,size(POS,3)),POS);
                    POS = cumsum(POS,1);
                    POS = squeeze(POS);
                    POS(:,end) = [];
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % get the velocity values
                    VEL = diff(Xwin,1,3);
                    VEL = sum(VEL.*VEL,2).^.5;
                    VEL = squeeze(VEL);
                    VEL = VEL * sam^-1;
                    
                    %{
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Xm
                    Xwin = tmp_Xm;
                    for loop = 1:100
                        for tm = 1:size(POS,2)
                            plot(POS(:),VEL(:),'b.')
                            hold on;
                            plot(POS(:,tm),VEL(:,tm),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);                            
                            hold off
                            drawnow
                            title(num2str(tm));
                            pause(.4)

                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %}
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % fit-interpolate data
                    tmpVEL = fitFun(POS(:),VEL(:),tmpDOMAIN,100,0);                    
                    tmpREGR = cwt(-tmpVEL,cwtSMOOTH,'gaus1');
                    dS = cwt(-tmpDOMAIN,cwtSMOOTH,'gaus1');
                    tmpREGR = tmpREGR.*dS.^-1;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % store the results
                    totalREGR = [totalREGR;tmpREGR];
                    totalVEL = [totalVEL;tmpVEL];
                    totalDomain1 = [totalDomain1;tmpDOMAIN];
                    totalDomain2 = [totalDomain2;seg*ones(size(tmpDOMAIN))];
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % plot results
                    if disp
                        plot(POS(:),VEL(:),'r.')
                        hold on;
                        plot(tmpDOMAIN,tmpVEL)
                        hold off
                        drawnow
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % PLOT: each as steady state
                    figure;                    
                    fidx = find(~isnan(tmpVEL));
                    % plot "raw"
                    plot(POS(:),VEL(:),'b.');
                    hold on
                    plot(tmpDOMAIN(fidx),tmpVEL(fidx),'r.');
                    title(['{[sample].[' num2str(sam) ']}{[window].[' num2str(WINDOW(win)) ']}{[seg].[' num2str(segValue) ']}']);
                    drawnow
                    saveas(gca,[outPath filesep '{[name].[' nm ']}{[sample].[' num2str(sam) '][window].[' num2str(win) ']}{[type][velocity-2D]}.tif']);
                    close all;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % store the results
                SNIP = 3*cwtSMOOTH;
                didx = find(~any(isnan(totalVEL),1));
                didx = didx(SNIP:end);                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % MESH: store the results - for non-steady state
                figure;mesh(totalDomain1(:,didx),totalDomain2(:,didx),totalVEL(:,didx));
                view([0 90]);
                saveas(gca,[outPath filesep '{[name].[' nm ']}{[sample].[' num2str(sam) '][window].[' num2str(win) ']}{[type][velocity]}.tif']);
                close all
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % MESH: store the results
                figure;mesh(totalDomain1(:,didx),totalDomain2(:,didx),totalREGR(:,didx));
                view([0 90]);
                saveas(gca,[outPath filesep '{[name].[' nm ']}{[sample].[' num2str(sam) '][window].[' num2str(win) ']}{[type][REGR]}.tif']);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % store the results
                csvwrite([outPath filesep '{[name].[' nm ']}{[sample].[' num2str(sam) '][window].[' num2str(win) ']}{[type][raw_REGR]}.csv'],totalREGR(:,didx));
                csvwrite([outPath filesep '{[name].[' nm ']}{[sample].[' num2str(sam) '][window].[' num2str(win) ']}{[type][raw_velocity]}.csv'],totalVEL(:,didx));
                close all
            end     
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store the results
        rawPoints = permute(PointList,[3 2 1]);
        rawPoints = reshape(rawPoints,[size(rawPoints,1) size(rawPoints,2)*size(rawPoints,3)]);
        csvwrite([outPath filesep '{[name].[' nm ']}{[sample].[' num2str(sam) '][window].[' num2str(win) ']}{[type][raw_rawpoints]}.csv'],rawPoints);
    
        
        %}
        
        waitforbuttonpress
    catch ME
        ME
        
        
    end
    
end
%%
e = 20;
SCALE = .25;
tmpI = double(I(:,:,e))/255;
N = round([100 500]*SCALE);
R = round([0 200]*SCALE);
tmpO = tmpI;
tmpI = imresize(tmpI,SCALE);
tmpR = tmpI;
for r = 1:4
    tmpR(:,1:R(2)) = [];
    tmpR = imrotate(tmpR,90);
end

[n1 n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
[d1 d2] = ndgrid(R(2):(size(tmpI,1)-R(2)),R(2):(size(tmpI,2)-R(2)));
X = n1.*cos(n2);
Y = n1.*sin(n2);
close all
disp = 0;
g1 = {};
g2 = {};
parfor p = 1:numel(d1)
    
    LOC = [d1(p) d2(p)];
    %[LOC(1) LOC(2) V] = impixel(tmpI);
    Xp = X + LOC(2);
    Yp = Y + LOC(1);
    if disp
        imshow(tmpI);
        hold on
        for s = 1:1:size(Xp,2)
            plot(Xp(:,s),Yp(:,s),'r');
        end
        for s = 1:1:size(Xp,1)
            plot(Xp(s,:),Yp(s,:),'r');
        end
        hold off
        drawnow
    end
    F = ba_interp2(tmpI,Xp,Yp);
    F = bsxfun(@minus,F,mean(F,2));
    fF = fft(F,[],2);
    ufF = mean(abs(fF),1);
    [q1 q2] = sort(ufF(1:end/2),'descend');
    g1{p} = q1(1:3);
    g2{p} = q2(1:3);
    %drawnow
    %hold off
    p
    numel(d1)
end
%%
g1 = cell2mat(g1');
g2 = cell2mat(g2');
%%
close all
T = reshape(g1(:,1:3),[size(d1) 3]);

imshow(T/max(T(:)),[]);
figure;imshow(tmpR);
L = g2(:,1:3);
[UQ L L2] = unique(L,'rows');
L2 = label2rgb(L2);

T = reshape(L2,[size(d1) 3]);

figure;imshow(T);

T = reshape(g2(:,1:3),[size(d1) 3]);
figure;imshow(T(:,:,1) == 3);
M = T(:,:,1) == 3;

M = imresize(M,size(tmpR));
Re = regionprops(M,'MajorAxisLength','PixelIdxList');
[J midx] = max([Re.MajorAxisLength]);
M = zeros(size(M));
M(Re(midx).PixelIdxList) = 1;
M = imfill(M,'holes');

M = padarray(M, [R(2) R(2)]);
M = imresize(M,size(tmpO));

M = imclose(logical(M),strel('disk',50));
skel = bwmorph(M,'thin',inf);
skel = imdilate(skel,strel('disk',3));

figure;
imshow(M,[]);
out = flattenMaskOverlay(tmpO,logical(skel));
figure;
imshow(out);