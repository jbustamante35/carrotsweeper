%% get list of data to generate "synthetic" hand data
dataPath = ['/iplant/home/hirsc213/maizeData/seedlingData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
FileList = {};
FileExt = {'nef'};
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
% submit N files for training networks
N = 500;
subList = FileList(1:N);
[subList] = issueBulkTicket(subList);
reSize = .25;
newY = 1200;
% spin up cFlow
func = cFlow('getTrainData');
func.setMCRversion('v840');
func.setMemory(8000);
numJobs = N;
% init vars for holding pointers
MASK = {};
sI = {};
sM = {};
rI = {};
% local run
% getTrainData(subList{15},1200);
for e = 1:N
    tic;
    [MASK{e},sI{e},sM{e},rI{e}] = func(subList{e},newY);
    toc
    e
end
% submit data to condor
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,500,500);
%% load images and mask for making vertical strip choice
I = imread(FileList{1});
[I angle] = rectifyImage(I);
I = imresize(I,.25);
mN = 500;
Wl = zeros([size(I) mN]);
MASKl = cell(mN,1);
for e = 1:mN
    try
        Wl(:,:,:,e) = cFlowLoader(rI{e});
        MASKl{e} = cFlowLoader(MASK{e});
        kp(e) = 1;
        e
    catch
        kp(e) = 0;
    end
end
% stack for vertical strips
Y = [];
for e = 1:mN
    try
        tmp = full(MASKl{e});
        tmp = imresize(tmp,reSize,'nearest');
        tmp = any(tmp,1);
        Y = [Y tmp];
        e
    catch
        kp(e) = 0;
    end
end
%%
Wl(:,:,:,find(kp==0)) = [];
%% train vertical strip network
[funcObject1] = getNetwork(Wl,5,15,Y);
funcObject1_functionName = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/smartA/maizeSeedling_func1.m';
genFunction(funcObject1.net,funcObject1_functionName);
funcObject1.func = @(X)maizeSeedling_func1(X);
clear Wl Y
%% stack for horizontal strips
Y2 = [];
X2 = [];
%Y2 = zeros([1 size(sM{1},1)*numel(sM)]);
toStack = 300;
%X2 = zeros([size(sIl{1},2) size(sIl{1},1) 3 3*toStack]);
cnt = 1;
parfor e = 1:toStack
    try 
        tmpM{e} = cFlowLoader(sM{e});
        tmpI{e} = cFlowLoader(sI{e});
         kp2(e) = 1;
    catch
        kp2(e) = 0;
    end
    e
end

%% permute and stack
cnt = 1;
X2 = zeros([size(tmpI{1},2) size(tmpI{1},1) 3 3*toStack]);
Y2 = [];
for e = 1:toStack
    for r = 1:3
        Y2 = [Y2 any(tmpM{e}(:,:,r),2)'];
        X2(:,:,:,cnt) = permute(tmpI{e}(:,:,:,r),[2 1 3]);
        cnt = cnt + 1;
    end
    e
end
%% try to conver to t hsv
for e = 1:size(X2,4)
    HX2(:,:,:,e) = rgb2hsv_fast(X2(:,:,:,e));
end
%%
[funcObject2] = getNetwork(sort(X2,1),5,10,Y2);
funcObject2_functionName = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/smartA/maizeSeedling_func2.m';
genFunction(funcObject2.net,funcObject2_functionName);
funcObject2.func = @(X)maizeSeedling_func2(X);
%%
cnt = 1;
%X2 = zeros([size(tmpI{1},2) size(tmpI{1},1) 3 3*toStack]);
%Y2 = [];
for e = 1:toStack
    for r = 1:3
        fidx = find(any(tmpM{e}(:,:,r),2)');
        subI = tmpI{e}(fidx,:,:,r);
        X3{e} = permute(subI,[2 1 3]);
        %imshow(X3{e},[])
        %drawnow
        cnt = cnt + 1;
    end
    e
end
%%
B = [];
for e = 1:numel(X3)
    B = [B X3{e}];
end
%%
B = sort(B,1);
%B = permute(B,[2 1 3 4]);
B = permute(B,[1 3 2 4]);
sz = size(B);
B = reshape(B,[prod(sz(1:2)) prod(sz(3))]);
[bS bC bU bE bL bERR bLAM] = PCA_FIT_FULL_T(B,5);

%%
obj = gmdistribution.fit(bC',3);
kidx = cluster(obj,bC');
%%
[funcObject3] = getNetwork(double(B),2,10,Y3);
funcObject3_functionName = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/smartA/maizeSeedling_func3.m';
genFunction(funcObject3.net,funcObject3_functionName);
funcObject3.func = @(X)maizeSeedling_func3(X);
%%
Y3 = [];
str = 1;
for e = 1:numel(X3)
    stp = str + size(X3{e},2) - 1;
    tmpM = kidx(str:stp)'==3;
    R = regionprops(tmpM,'PixelIdxList');
    tmpM = zeros(size(tmpM));
    tmpM(R(end).PixelIdxList) = 1;
    tmpM = logical(bwlarge(tmpM));
    Y3 = [Y3 tmpM];
    str = stp + 1;
    %plot(tmpM)
    %drawnow
    %{
    out = flattenMaskOverlay(X3{e},repmat(tmpM,[size(X3{e},1) 1]));
    imshow(out,[]);
    %}
    e
end



%{
%% 
Wl = zeros([size(I) N]);
MASKl = cell(1,N);
sIl = cell(1,N);
sMl = cell(1,N);
%%
for e = 650:N
    try
        Wl(:,:,:,e) = cFlowLoader(rI{e});
        MASKl{e} = cFlowLoader(MASK{e});
        sIl{e} = cFlowLoader(sI{e});
        sMl{e} = cFlowLoader(sM{e});
        kp(e) = 1;
        e
    catch
        kp(e) = 0;
    end
end
%}


%%
parfor e = 1:600
    BOT{e} = smartMain(FileList{e},reSize,funcObject1,funcObject2,funcObject3);
    e
end
%%
B = [];
for e = 1:numel(BOT)
    B = cat(4,B,BOT{e});
    e
    
    
    
    
    
end
B = sort(B,2);
B = permute(B,[2 1 3 4]);
B = permute(B,[1 3 2 4]);
sz = size(B);
B = reshape(B,[prod(sz(1:2)) prod(sz(3:4))]);
[bS bC bU bE bL bERR bLAM] = PCA_FIT_FULL_T(B,5);
%%
obj = gmdistribution.fit(bC',3);
kidx = cluster(obj,bC');
%%
close all
plot(kidx(1:200))
%%
X3 = [];
Y3 = [];
for e = 1:200
    tmpO = smartMain(FileList{e},reSize,funcObject1,funcObject2);

    tmp = sort(tmpO,2);
    tmp = permute(tmp,[2 1 3 4]);
    tmp = permute(tmp,[1 3 2 4]);
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:2)) prod(sz(3:4))]);
    tmp = PCA_REPROJ_T(tmp,bE,bU);
    tmp = reshape(tmp,[size(tmp,1) size(tmpO,1) size(tmpO,4)]);
    for i = 1:size(tmp,3)
        
        Yp = cluster(obj,tmp(:,:,i)');
        Y3 = cat(2,Y3,Yp'==3);
        sig = repmat((Yp==3),[1 size(tmpO,2)]);
        out = flattenMaskOverlay(tmpO(:,:,:,i)/255,sig);
        imshow(out,[]);
        drawnow
    end
    X3 = cat(4,X3,tmpO);
    
end
X3 = permute(X3,[2 1 3 4]);
%%
sum(kp==0)
%%
sz = size(W);
rm = reshape(W,[prod(sz(1:3)) sz(4)]);
rm = find(sum(rm,2)==0);


%%
Y = [];
for e = 1:numel(MASK)
    tmp = full(MASK{e});
    tmp = imresize(tmp,reSize,'nearest');
    tmp = any(tmp,1);
    Y = [Y tmp];
    e
end
%%%
sz = size(W);
W = permute(W,[1 3 2 4]);
sz1 = size(W);
W = reshape(W,[prod(sz1(1:2)) prod(sz1(3:4))]);
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL_T(W,5);
%%
net = patternnet(10);
net = train(net,wC,logical(Y),'useParallel','yes');
%%
wCr = reshape(wC,[5 sz1(3) sz1(4)]);
for e = 1:size(wCr,3)
    sig = net(wCr(:,:,e));
    sig = logical(sig > .5);
    sig = imclose(sig,ones([1 50]));
    R = regionprops(logical(sig),'Area','PixelIdxList');
    R([R.Area] < 100) = [];
    sig = zeros(size(sig));
    for r = 1:numel(R)
        sig(R(r).PixelIdxList) = 1;
    end
    
    I = imread(FileList{e});
    sig = imresize(sig,[1 size(I,2)]);
    sig = repmat(sig,[size(I,1) 1]);
    tmp = full(MASK{e});
    out = flattenMaskOverlay(double(I)/255,logical(tmp),.5,'r');
    out = flattenMaskOverlay(out,logical(sig>.5),.25,'b');
end
%% stack the other way
parfor e = 1:500
    
end
%% quick scan
for e = 1:numel(FileList{e})
    I = imread(FileList{e});
    [I angle] = rectifyImage(I);
    out = flattenMaskOverlay(I,MASK{e});
    imshow(out,[]);
    drawnow
end
%%
N = 500;
I = cell(N,1);
eBLOCK = cell(N,1);
parfor e = 1:N
    try
        I{e} = imread(FileList{e});
        [I{e} angle] = rectifyImage(I{e});
        eBLOCK{e} = findVerticalStrips(I{e});
        e
    catch
    end
end
%% crop vertical strip
for e = 1:numel(eBLOCK)
    [R{e}] = getMaskFromVerticalStrip(eBLOCK{e},size(I{e}),50,20);
    [tmpD{e}] = cropAndtrim(I{e},R{e},100);
    e
end

%%
ms = [];
msk = [];
H = 3325 - 1000;
for e = 1:numel(tmpD)
    for r = 1:numel(tmpD{e})
        ms = cat(4,ms,imresize(tmpD{e}{r}((end-H):end,:,:),[size(tmpD{e}{r}((end-H):end,:,:),1) 1200]));
        %msk= cat(3,msk,getMASK_ver0(imresize(tmpD{e}{r}((end-H):end,:,:),[size(tmpD{e}{r}((end-H):end,:,:),1) 1200])));
    end
    e
end
ms = permute(ms,[2 1 3 4]);
ms = double(ms)/255;
%% 
for e = 1:size(ms,4)
    ms(:,:,:,e) = imfilter(ms(:,:,:,e),fspecial('gaussian',[31 31],7),'replicate');
end
%% for hori vects

ms = permute(ms,[1 2 4 3]);
ms = sort(ms,1);
%dms = gradient(ms);
sz = size(ms);
ms = reshape(ms,[sz(1) prod(sz(2:3)) sz(4)]);
uM = mean(ms,2);
ms = bsxfun(@minus,ms,uM);
%%
sz2 = size(ms);
ms = reshape(ms,[sz(1) prod(sz2(2:3))]);
%%
[wS2 wC2 wU2 wE2 wL2 wERR2 wLAM2] = PCA_FIT_FULL_T(ms,5);
%%
pC = reshape(wC2,[5 sz2(2) sz2(3)]);
pC = reshape(pC,[5 sz(2:4)]);
pC = ipermute(pC,[1 2 4 3]);

%%
pC = permute(pC,[1 3 2 4]);
pC = reshape(pC,[size(pC,1)*size(pC,2) size(pC,3) size(pC,4)]);
pC = reshape(pC,[size(pC,1) size(pC,2)*size(pC,3)]);
%%

%%
close all
H = 3325 - 1000;
siG = [];
for e = 1:numel(tmpD)
    for r = 1:numel(tmpD{e})
         sig = zeros(2326,1);
        try
            Hi = double(imresize(tmpD{e}{r}((end-H):end,:,:),[size(tmpD{e}{r}((end-H):end,:,:),1) 1200]))/255;
            Hi = imfilter(Hi,fspecial('gaussian',[31 31],7),'replicate');
            E = edge(rgb2gray(Hi));
            E = imdilate(E,strel('disk',5,0));
            % find the 90-plumb
            [Hu,T,R] = hough(E','Theta',linspace(-5,5,100));
            P = houghpeaks(Hu,3);
            lines = houghlines(E',T,R,P,'FillGap',size(E,1)/2,'MinLength',500);
            xy = [lines(1).point1; lines(1).point2];
            v = round(mean(xy(:,1)));
            sig = zeros(size(Hi,1),1);
            sig(1:v) = 1;
            %Hi = double(tmpD{e}{r}((end-H):end,:,:))/255;
             Hi = double(imresize(tmpD{e}{r}((end-H):end,:,:),[size(tmpD{e}{r}((end-H):end,:,:),1) 1200]))/255;
            imshow(Hi,[]);
            hold on;
            plot(100*sig,1:numel(sig),'r')
            hold off
            drawnow
            siG = [siG sig];
        catch
            
        end
    end
end
%%
[wS3 wC3 wU3 wE3 wL3 wERR3 wLAM3] = PCA_FIT_FULL_T(pC,4);
GMModel = fitgmdist(wC3',4);
kidx = cluster(GMModel,wC3');
%%
close all
%kidx = kmeans(pC',3);
for e = 1:numel(tmpD)
    for r =1:3
      
        SKIP = size(tmpD{e}{r}((end-H):end,:,:),1);
        imshow(tmpD{e}{r}((end-H):end,:,:),[])
        hold on
        plot(100*kidx(1:1+SKIP-1),1:SKIP)
        hold off
        drawnow
    end
end
%%
IS = [];
IB = [];
cnt = 1;
for e = 1:numel(eBLOCK)
    if ~isempty(I{e})
        IS(:,:,:,cnt) = (double(I{e})/255);
        IB = [IB;(eBLOCK{e})];
        cnt = cnt + 1;
    end
    e
end
%% for vertical vects
IS = permute(IS,[1 2 4 3]);
sz = size(IS);
IS = reshape(IS,[sz(1) prod(sz(2:3)) sz(4)]);
uS = mean(IS,2);
IS = bsxfun(@minus,IS,uS);
%%
sz2 = size(IS);
IS = reshape(IS,[sz(1) prod(sz2(2:3))]);
%%
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL_T(IS,10);
%%
pC = reshape(wC,[10 sz2(2) sz2(3)]);
pC = reshape(pC,[10 sz(2:4)]);
pC = ipermute(pC,[1 2 4 3]);

%%
pC = permute(pC,[1 3 2 4]);
pC = reshape(pC,[size(pC,1)*size(pC,2) size(pC,3) size(pC,4)]);
pC = reshape(pC,[size(pC,1) size(pC,2)*size(pC,3)]);
%%
Yw = [zeros(163,30) IB zeros(163,71)];
[XL,YL,XS,YS,BETA,PCTVAR] = plsregress(pC',Y,15);
Yp = [ones(size(pC,2),1) pC']*BETA;
Yp = reshape(Yp,size(Yw'));
Yp = Yp';
%%
close all
for e = 1:size(Yp,1)
    plot(Yp(e,:));
    hold all
    plot(Yw(e,:)+.1);
    plot(Yp(e,:) > .4,'r')
    drawnow
    hold off
    pause(.1)
end
%%
ksdensity(Yp(find(Yw==1)))
%%
I = double(imread(FileList{3000}))/255;
I = bsxfun(@minus,I,uS);
pI = reshape(I,[size(I,1) size(I,2)*size(I,3)]);
[pC] = PCA_REPROJ_T(pI,wE,wU);
pC = reshape(pC,[10 size(I,2) size(I,3)]);
pC = permute(pC,[1 3 2]);
pC = reshape(pC,[size(pC,1)*size(pC,2) size(pC,3)]);




%%
e = 1;
%toM = [];
parfor e = 1:100
    try
        tic
        [mB{e}] = smartMain(FileList{e},.25,funcObject1,funcObject2,funcObject3,'./forCory/',[],toM);
        toc
        worked(e) = 1;
    catch ME
        worked(e) = 0;
    end
end
%%  publish NN

%%
toM = mB{1};
for e = 1:numel(mB)
    if worked(e)
        for m = 1:numel(mB{e})
            toM(m).statBlock = toM(m).statBlock + mB{e}(m).statBlock;
        end
    end
end

for m = 1:numel(toM)
    toM(m).statBlock = toM(m).statBlock/sum(toM(m).statBlock(:));
    toM(m).statBlock = -log(toM(m).statBlock);
end




