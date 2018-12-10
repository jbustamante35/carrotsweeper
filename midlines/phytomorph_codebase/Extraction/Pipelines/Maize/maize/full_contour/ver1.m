FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% remove non 61
ridx = [];
for e = 1:numel(SET)
    if numel(SET{e}) ~= 61
        ridx(e) = 1;
    else
        ridx(e) = 0;
    end
end
SET(find(ridx)) = [];
%% for each set
sampleBank = [];
curveBank = [];
ID = [];
cnt = 1;
for e = 1:5%numel(SET)
    tic
    %isolateKernels(SET{e},1,'/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction2/');
    data = sampleContours(SET{e},1,[]);    
    
    
    for e = 1:numel(data)
        hold off
        plot(data(e).contourdata.data(1,:),data(e).contourdata.data(2,:),'k')
        hold on
        plot(data(e).radical.data(1,:),data(e).radical.data(2,:),'r.')
        hold off
        pause(.1)
        
    end
    
    
    toc
end
%% group and view banks
[S C U E L ERR LAM] = PCA_FIT_FULL(sampleBank',10);
NG = 4;
kidx = kmeans(C,NG);
UQ = unique(ID);
CL = {'r.' 'b.'};
for u = 1:numel(UQ)
    fidx = find(ID == UQ(u));
    tc = curveBank(:,fidx);
    tk = kidx(fidx);
    plot(tc(1,:),tc(2,:),'k')
    hold on
    for i = 1:NG
        plot(tc(1,find(tk==i)),tc(2,find(tk==i)),CL{i});
    end    
    drawnow
    
end
%%
I = imread(SET{1}{1});
S = zeros([size(I) numel(SET)]);
e
rm = [];
for e = 1:numel(SET)
        try
        S(:,:,e) = imread(SET{e}{1});
        %imshow(S(:,:,e),[])
        %drawnow
        e
        catch
            rm = [rm e];
        end
end
%%
S(:,:,rm) = [];
%%
% for non max
fS2 = 31;
fS1 = 201;
fS = 31;
gradDiskSize = 11;
% for plots
mag2 = 500;
mag1 = 3000;
% row samples width
rowS = 10;
curveBank = {};
for e = 1:size(S,3)
    try
        I = S(10:end-10,300:end-500,e);
        BK = imclose(double(I),strel('disk',51));
        I = I - BK;
        I = bindVec(I);
        
        
        [d1 d2] = gradient(imfilter(I,fspecial('disk',gradDiskSize),'replicate'));
        G = (d1.*d1 + d2.*d2).^.5;
        G = bindVec(G);
        thresholdG = graythresh(G);
        E = G > thresholdG;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prepare for rows
        Ii = abs(I-1);
        
        % find rows
        %s2 = std(abs(d1),1,2);
        %s2 = sum(I,2);
        s2 = sum(Ii,2).*std(abs(d1),1,2);
        s2 = bindVec(s2);
        s2 = imfilter(s2,fspecial('disk',fS),'replicate');
        %es2 = imerode(s2,ones(fS2,1));
        es2 = imdilate(s2,ones(fS2,1));
        p2 = s2 == es2;

        fidx = find(p2);
        sam2 = s2(fidx);
        thresh2 = graythresh(s2);
        %fidx = fidx(sam2 < thresh2);
        fidx = fidx(sam2 > thresh2);
        p2 = zeros(size(p2));
        p2(fidx) = 1;
        P2 = repmat(p2,[1 size(I,2)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find cols
        s1 = sum(Ii,1).*std(abs(d2),1,1);
        %s1 = sum(I,1);
        s1 = bindVec(s1);
        s1 = imfilter(s1,fspecial('disk',fS),'replicate');
        %es1 = imerode(s1,ones(1,fS1));
        es1 = imdilate(s1,ones(1,fS1));
        p1 = s1 == es1;

        fidx = find(p1);
        sam1 = s1(fidx);
        thresh1 = graythresh(s1);
        %fidx = fidx(sam1 < thresh1);
        %fidx = fidx(sam1 > thresh1);
        [J sidx] = sort(sam1);
        p1 = zeros(size(p1));
        p1(fidx(sidx(end))) = 1;
        P1 = repmat(p1,[size(I,1) 1]);

      
        
        
        %%% find plug sites
        %fM = levelSetKurvatureFeatureMap();
        %fM = fM.computeFeatureMap(imfilter(I,fspecial('disk',11),'replicate'));
        
        CENTERS = P1.*P2;
        
        % FIND THE CAP
        % sample along each row
        clear cI CROPBOX;
        fidx2 = find(p2);
        for r = 1:sum(p2)
            rowSample = E(fidx2(r)-rowS:fidx2(r)+rowS,:);
            rowSample = mean(rowSample,1);
            capIdx = find(rowSample ~= 0);
            cI{r} = [capIdx(1) + (gradDiskSize-1)/2   fidx2(r)];
            UL = cI{r} - [20 100];
            BR = cI{r} + [500 100];
            SZ = BR - UL;
            CROPBOX{r} = [UL SZ];
        end
        [cy cx] = find(CENTERS);
        
        
        % FIND THE TOP AND BOTTOM
        % sample at each center
        TOP = {};
        BOTTOM = {};
        for r = 1:numel(cx)
            colSampleUp = E(1:cy(r),cx(r)-10:cx(r)-10);
            colSampleUp = mean(colSampleUp,2);
            colSampleDown = E(cy(r):end,cx(r)-10:cx(r)-10);
            colSampleDown = mean(colSampleDown,2);
            topIDX = find(colSampleUp);
            bottomIDX = find(colSampleDown);
            TOP{r} = [cx(r) topIDX(end) - (gradDiskSize-1)/2];
            BOTTOM{r} = [cx(r) cy(r)+bottomIDX(1) + (gradDiskSize-1)/2];
        end
        
        
        imshow(I,[])
        hold on;
        plot(mag2*p2,1:size(I,1),'r');
        plot(1:size(I,2),mag1*p1,'g');
        %
        for i = 1:numel(cI)
            plot(cI{i}(1),cI{i}(2),'b*')
        end
        %
        for i = 1:numel(TOP)
            plot(TOP{i}(1),TOP{i}(2),'c*')
            plot(BOTTOM{i}(1),BOTTOM{i}(2),'y*')
        end
        %
        for i = 1:numel(TOP)
            plot(cx(i),cy(i),'m*')
        end
        % analysis of each image box
        for i = 1:numel(CROPBOX)
            sI{i} = imcrop(I,CROPBOX{i});
            sI{i} = imfilter(sI{i},fspecial('gaussian',11,5),'replicate');
            hold on;
            C = contourc(sI{i},10);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % get the curve structure
            str = 1;
            c = 1;
            clear curve
            while str < size(C,2)
                ed = str + C(2,str);
                curve(c).level = C(1,str);
                curve(c).data = C(:,str+1:ed);
                curve(c).length = size(curve(c).data,2);
                c = c + 1;
                str = ed + 1;
            end
            [curve] = selectClosedCurves(curve);
            imshow(sI{i})
            [J,sidx] = sort([curve.length]);
            for j = 1:numel(curve)
                plot(curve(j).data(1,:),curve(j).data(2,:),'r');
            end
            plot(curve(sidx(end)).data(1,:),curve(sidx(end)).data(2,:),'g','LineWidth',2);
            out = cwtK_closed(curve(sidx(end)).data',{60});
            [J tipIDX] = min(out.K);
            plot(curve(sidx(end)).data(1,tipIDX),curve(sidx(end)).data(2,tipIDX),'k*');
            % shift for tipcentered
            cur = circshift(curve(sidx(end)).data(:,1:end-1),[0 -tipIDX]);
            plot(cur(1,1),cur(2,1),'ko');
            %out = cwtK_closed(cur',{5});
            curveBank{end+1} = cur;
            drawnow
            pause(.1);
% for non max
fS2 = 31;
fS1 = 201;
fS = 31;
gradDiskSize = 11;
% for plots
mag2 = 500;
mag1 = 3000;
% row samples width
rowS = 10;
curveBank = {};
for e = 1:size(S,3)
    try
        I = S(10:end-10,300:end-500,e);
        BK = imclose(double(I),strel('disk',51));
        I = I - BK;
        I = bindVec(I);
        
        
        [d1 d2] = gradient(imfilter(I,fspecial('disk',gradDiskSize),'replicate'));
        G = (d1.*d1 + d2.*d2).^.5;
        G = bindVec(G);
        thresholdG = graythresh(G);
        E = G > thresholdG;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prepare for rows
        Ii = abs(I-1);
        
        % find rows
        %s2 = std(abs(d1),1,2);
        %s2 = sum(I,2);
        s2 = sum(Ii,2).*std(abs(d1),1,2);
        s2 = bindVec(s2);
        s2 = imfilter(s2,fspecial('disk',fS),'replicate');
        %es2 = imerode(s2,ones(fS2,1));
        es2 = imdilate(s2,ones(fS2,1));
        p2 = s2 == es2;

        fidx = find(p2);
        sam2 = s2(fidx);
        thresh2 = graythresh(s2);
        %fidx = fidx(sam2 < thresh2);
        fidx = fidx(sam2 > thresh2);
        p2 = zeros(size(p2));
        p2(fidx) = 1;
        P2 = repmat(p2,[1 size(I,2)]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find cols
        s1 = sum(Ii,1).*std(abs(d2),1,1);
        %s1 = sum(I,1);
        s1 = bindVec(s1);
        s1 = imfilter(s1,fspecial('disk',fS),'replicate');
        %es1 = imerode(s1,ones(1,fS1));
        es1 = imdilate(s1,ones(1,fS1));
        p1 = s1 == es1;

        fidx = find(p1);
        sam1 = s1(fidx);
        thresh1 = graythresh(s1);
        %fidx = fidx(sam1 < thresh1);
        %fidx = fidx(sam1 > thresh1);
        [J sidx] = sort(sam1);
        p1 = zeros(size(p1));
        p1(fidx(sidx(end))) = 1;
        P1 = repmat(p1,[size(I,1) 1]);

      
        
        
        %%% find plug sites
        %fM = levelSetKurvatureFeatureMap();
        %fM = fM.computeFeatureMap(imfilter(I,fspecial('disk',11),'replicate'));
        
        CENTERS = P1.*P2;
        
        % FIND THE CAP
        % sample along each row
        clear cI CROPBOX;
        fidx2 = find(p2);
        for r = 1:sum(p2)
            rowSample = E(fidx2(r)-rowS:fidx2(r)+rowS,:);
            rowSample = mean(rowSample,1);
            capIdx = find(rowSample ~= 0);
            cI{r} = [capIdx(1) + (gradDiskSize-1)/2   fidx2(r)];
            UL = cI{r} - [20 100];
            BR = cI{r} + [500 100];
            SZ = BR - UL;
            CROPBOX{r} = [UL SZ];
        end
        [cy cx] = find(CENTERS);
        
        
        % FIND THE TOP AND BOTTOM
        % sample at each center
        TOP = {};
        BOTTOM = {};
        for r = 1:numel(cx)
            colSampleUp = E(1:cy(r),cx(r)-10:cx(r)-10);
            colSampleUp = mean(colSampleUp,2);
            colSampleDown = E(cy(r):end,cx(r)-10:cx(r)-10);
            colSampleDown = mean(colSampleDown,2);
            topIDX = find(colSampleUp);
            bottomIDX = find(colSampleDown);
            TOP{r} = [cx(r) topIDX(end) - (gradDiskSize-1)/2];
            BOTTOM{r} = [cx(r) cy(r)+bottomIDX(1) + (gradDiskSize-1)/2];
        end
        
        
        imshow(I,[])
        hold on;
        plot(mag2*p2,1:size(I,1),'r');
        plot(1:size(I,2),mag1*p1,'g');
        %
        for i = 1:numel(cI)
            plot(cI{i}(1),cI{i}(2),'b*')
        end
        %
        for i = 1:numel(TOP)
            plot(TOP{i}(1),TOP{i}(2),'c*')
            plot(BOTTOM{i}(1),BOTTOM{i}(2),'y*')
        end
        %
        for i = 1:numel(TOP)
            plot(cx(i),cy(i),'m*')
        end
        % analysis of each image box
        for i = 1:numel(CROPBOX)
            sI{i} = imcrop(I,CROPBOX{i});
            sI{i} = imfilter(sI{i},fspecial('gaussian',11,5),'replicate');
            hold on;
            C = contourc(sI{i},10);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % get the curve structure
            str = 1;
            c = 1;
            clear curve
            while str < size(C,2)
                ed = str + C(2,str);
                curve(c).level = C(1,str);
                curve(c).data = C(:,str+1:ed);
                curve(c).length = size(curve(c).data,2);
                c = c + 1;
                str = ed + 1;
            end
            [curve] = selectClosedCurves(curve);
            imshow(sI{i})
            [J,sidx] = sort([curve.length]);
            for j = 1:numel(curve)
                plot(curve(j).data(1,:),curve(j).data(2,:),'r');
            end
            plot(curve(sidx(end)).data(1,:),curve(sidx(end)).data(2,:),'g','LineWidth',2);
            out = cwtK_closed(curve(sidx(end)).data',{60});
            [J tipIDX] = min(out.K);
            plot(curve(sidx(end)).data(1,tipIDX),curve(sidx(end)).data(2,tipIDX),'k*');
            % shift for tipcentered
            cur = circshift(curve(sidx(end)).data(:,1:end-1),[0 -tipIDX]);
            plot(cur(1,1),cur(2,1),'ko');
            %out = cwtK_closed(cur',{5});
            curveBank{end+1} = cur;
            drawnow
            pause(.1);
        end
         
        
        drawnow
        hold off
    catch ME
        fprintf(['error@' num2str(e) '\n']);
    end
end
        end
         
        
        drawnow
        hold off
    catch ME
        fprintf(['error@' num2str(e) '\n']);
    end
end
%%
N = 301;
CK = [];
for e = 1:numel(curveBank)
    
    tmpCurve = curveBank{e};
    tmpCurve = bsxfun(@minus,tmpCurve,tmpCurve(:,1));
    dL = diff(tmpCurve,1,2);
    L = cumsum([0 sum(dL.*dL,1).^.5]);
    H = interp1(L',tmpCurve',linspace(0,L(end),N)');
    
    
    if e == 1
        meanCurve = H;
    else
        meanCurve = .5*(meanCurve + H);
        %meanCurve = .5*(meanCurve + [H(:,1);-H(:,2)]);
    end
    
    fn{e} = spap2(20,5,L,tmpCurve);
    LEN(e) = L(end);
    
    plot(tmpCurve(1,:),tmpCurve(2,:),'r','LineWidth',2);
    hold on;
    fdata = fnval(fn{e},L);
    
    %plot(tmpCurve(1,:),tmpCurve(2,:),'k');
    %hold off
    %drawnow
    %pause(.3)
    out = cwtK_closed(H,{5});
    CK = [CK;out.K];
end
%%
close all
dL = diff(meanCurve,1,1);
L = cumsum([0;sum(dL.*dL,2).^.5]);
ufn = spap2(20,5,L',meanCurve');
UL = L(end);
curBase = fnval(ufn,linspace(0,UL,UL));
plot(curBase(1,:),curBase(2,:));
hold on;
axis equal
plot(curBase(1,round(UL/2)),curBase(2,round(UL/2)),'r*');
plot(curBase(1,round(UL/4)),curBase(2,round(UL/4)),'r*');
plot(curBase(1,round(UL*3/4)),curBase(2,round(UL*3/4)),'r*');
%% 
SAMREG_number = N-1;
POP = 700;
GENS = 1000;
MS = 20;
% set the optimization parameters
options = psooptimset('ConstrBoundary','reflect','PopInitRange',[zeros(1,SAMREG_number);MS*ones(1,SAMREG_number)],'PopulationSize',POP,'Display','on','Generations',GENS,'Vectorized','on','CognitiveAttraction',1.5,'SocialAttraction',.5);    
% for tangnet alignment
tr = 20;
% run the optimization
%xo = pso(@(xo)curveDistance2(ufn,fn{tr},UL,LEN(tr),xo),SAMREG_number,[],[],[],[],zeros(1,SAMREG_number),MS*ones(1,SAMREG_number),[],options);

% for kurvature alignment
% run the optimization
xo = pso(@(xo)curveDistance(fn{tr},LEN(tr),xo,cumsum(mean(CK,1))),SAMREG_number,[],[],[],[],zeros(1,SAMREG_number),MS*ones(1,SAMREG_number),[],options);

%
L = [0 cumsum(xo,2)];
L = bsxfun(@times,L,L(:,end).^-1);
L = L*LEN(tr);

CURV = fnval(fn{tr},L);
CURVi = fnval(fn{tr},linspace(0,LEN(tr),LEN(tr)));

per = .15;
hold on;
plot(CURVi(1,:),CURVi(2,:),'k');

plot(CURV(1,:),CURV(2,:),'r*');
plot(CURV(1,.5*SAMREG_number),CURV(2,.5*SAMREG_number),'g*');
plot(CURV(1,per*SAMREG_number),CURV(2,per*SAMREG_number),'k*');
plot(CURV(1,(1-per)*SAMREG_number),CURV(2,(1-per)*SAMREG_number),'k*');
plot(CURV(1,1),CURV(2,1),'g*');
axis equal

%%

N = 301;
CK = [];
for e = 1:numel(curveBank)
    
    tmpCurve = curveBank{e};
    tmpCurve = bsxfun(@minus,tmpCurve,tmpCurve(:,1));
    dL = diff(tmpCurve,1,2);
    L = cumsum([0 sum(dL.*dL,1).^.5]);
    H = interp1(L',tmpCurve',linspace(0,L(end),N)');
    
    cp = mean(H,1);
end
matched_lower = boundaryManifold(mag_lower,mag_upper,GEN,POP,0);

