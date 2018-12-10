%% look for many movies
FilePath = '/mnt/tetra/nate/confocalAlign/';
FileList = {};
FileExt = {'lsm'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% look for many movies
FilePath = '/mnt/tetra/nate/Nate/';
FileList = {};
FileExt = {'czi'};
FileList = gdig(FilePath,FileList,FileExt,1);

%% view movies
close all
h1= figure;
h2 = figure;
sel = 2;

TAGNUM = numel(info)-1;
mag = 500;
disp = 0;
toP = [];

clear toP K
for sel = 1:numel(FileList)
    [toP{sel} K{sel}] = confocalRatio(FileList{sel},0,0,1);
end

%% view movies
close all
h1= figure;
h2 = figure;
sel = 2;

TAGNUM = numel(info)-1;
mag = 500;
disp = 0;
toP = [];


for sel = 1:numel(FileList)
    %info = imfinfo(FileList{sel});
    %N = (numel(info)/2);
    %toP = zeros(1,numel(info)/2);
    
    fs = [];
    data = bfopen(FileList{sel});
    for e = 1:(size(data{1},1)/3)
        idxB = (e-1)*3;
        for k = 1:3
            fs(:,:,k,e) = fliplr(data{1}{idxB+k,1});
        end
    end
    fs = double(fs);
    toP = [];
    K = [];
    parfor e = 1:size(fs,4)-5
        tic
        
        idx = (e-1)*2+1;
        %sourceO = double(imread(FileList{sel},[],idx))/255;
        sourceO = fs(:,:,:,e)/255;
        source = mean(sourceO(:,:,1:2),3);

        MASK = source > .6*graythresh(souce);



        squareMASK = zeros(size(source));


        sig1 = mean(source,1);
        sig1 = imfilter(sig1,fspecial('disk',11),'replicate');
        sig2 = imfilter(sig1,fspecial('disk',51),'replicate');
        gsig = gradient(sig2);
        [~,gidx] = max(sig2);
        [~,gidx2] = max(gsig);
        uT = mean(sig1((end-30):end));
        uT = uT*1.2;
        th1 = sig1 < uT;
        th1(1:gidx) = 0;
        th2 = zeros(size(th1));
        th3 = zeros(size(th1));
        th2(1:gidx) = 1;
        th3(1:gidx2) = 1;
        fidx = find(th1);
        th1 = zeros(size(th1));
        th1(1:fidx(1)) = 1;


        squareMASK = zeros(size(source));
        squareMASK(:,gidx2:fidx(1)) = 1;
        MASK = MASK.*squareMASK;
        out = flattenMaskOverlay(source,logical(MASK));

        sourceO = imfilter(sourceO,fspecial('disk',5),'replicate');
        RATIO = sourceO(:,:,2).*sourceO(:,:,1).^-1;
        toM = MASK.*RATIO;
        Y = sum(toM,1);
        X = sum(MASK,1);
        Y = Y.*X.^-1;
        xidx = find(X);
        Y = Y(xidx);
        k = interp1(Y,linspace(1,numel(Y),200));
        K(:,e) = k';
        toP(e) = sum(toM(:))/sum(MASK(:));
        %{
        if disp
            figure(h1);
            imshow(out,[])
            hold on
            plot(sig2*mag*10,'y')
            plot(th1*mag/10,'r')
            plot(th2*mag/10,'b');
            plot(th3*mag/10,'c');

            hold off
            figure(h2);
            plot(toP);
            drawnow
        end
        %}

        toc*N/10
    end
    SIG{sel} = toP;
    SIGK{sel} = K;
    figure(h1);
    plot(SIG{sel})
    hold all
    drawnow
end
%%

%%
close all
figure;
plot(toP);
%% remove any unwanted
FileList(4)= [];
%% self reg movies
regMovie = {};
regForm = {};
for e = 1:numel(FileList)
    fprintf(['Start with movie:' num2str(e) ':' num2str(numel(FileList)) '\n']);
    [regMovie{e} other{e} regForm{e}] = registerMovieToItself(FileList{e});
    fprintf(['Done with movie:' num2str(e) ':' num2str(numel(FileList)) '\n']);
end
%% remove any unwanted
regMovie(4) = [];
%% show movie means
close all
for e = 1:numel(regMovie)
    U = mean(regMovie{e},3);
    imshow(U,[]);
    drawnow
end
%% bundle for matlab folks
U1 = mean(regMovie{47},3);
for e = 1:numel(regMovie)
    U2{e} = mean(regMovie{e},3);
end
%% RUN THIS
close all
regImage = {};
regFormBetween = {};
disp = 1;
for e = 3:numel(U2)
    if disp
        %RGB = cat(3,U1,U2{e},zeros(size(U1)));
        %imshow(RGB,[]);
        %drawnow
    end
    [regImage{e} regFormBetween{e} ltForm{e}] = regBetweenMovie(U1,U2{e},0);
    if disp
        %RGB = cat(3,U1,mean(U2{e},3),zeros(size(U1)));
        RGB = regImage{e};
        imshow(RGB,[]);
        title('done');
        drawnow
        pause(2);
        %waitforbuttonpress
    end
end
%% reg movies
close all
regImage = {};
regFormBetween = {};
disp = 1;
for e = 1:numel(regMovie)
    U1 = mean(regMovie{47},3);
    U2 = mean(regMovie{e},3);
    %{
    if disp
        RGB = cat(3,U1,U2,zeros(size(U1)));
        imshow(RGB,[]);
        drawnow
    end
    %}
    [regImage{e} regFormBetween{e} ltForm{e}] = regBetweenMovie(U1,U2,0);
    if disp
        %U2 = mean(regImage{e},3);
        %RGB = cat(3,U1,U2,zeros(size(U1)));
        RGB = regImage{e};
        imshow(RGB,[]);
        title('done');close
        drawnow
        pause(2);
        %waitforbuttonpress
    end
end
%% apply all transformations
udata = [0 1];  vdata = [0 1];
Rout = imref2d(size(mean(regMovie{1},3)));
Rout.XWorldLimits(2) = Rout.XWorldLimits(2);
Rout.YWorldLimits(2) = Rout.YWorldLimits(2);
Rout.ImageSize = Rout.ImageSize;
for e = 1:numel(FileList)
    rStack{e} = grandWarp(FileList{e},regForm{e},regFormBetween{e},ltForm{e},Rout,sM);
    fprintf(['Done with stack:' num2str(e) ':' num2str(numel(FileList)) '\n']);
end
%% stack registered into square frame
for e = 1:numel(rStack)
    SZ(e) = size(rStack{e},4);
end
ms = min(SZ);
mySZ = size(rStack{1});
mySZ(4) = ms;
MOV = zeros(mySZ);
for e = 1:numel(rStack)
    e
    MOV(:,:,:,:,e) = rStack{e}(:,:,:,1:ms);
end
%% remove blue channel from MOV
MOV(:,:,3,:,:) = [];
%% reshape for decompose MOV
SZ = size(MOV);
tot = reshape(MOV,[SZ(1:3) prod(SZ(4:5))]);
SZ1 = size(tot);
tot = reshape(tot,[prod(SZ1(1:3)) SZ1(end)]);
size(tot)
%% decompose
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL_T(tot,3);
%% reshape sim movies
SIM = reshape(tS,[SZ1]);
SIM = reshape(SIM,SZ);
%% reshape coeffs
ttC = reshape(tC,[3 SZ1(4)]);
%% reshape coeffs
ttC = reshape(ttC,[3 SZ(4) SZ(5)]);
%% 
sU = mean(sMOV,5);
%%
R = MOV(:,:,1,:,:).*MOV(:,:,2,:,:).^-1;
%%
U1 = squeeze(mean(R,5));
%%  
MASK = mean(mean(MOV(:,:,1,:,:),5),4);
%%
close all
M = MASK > graythresh(MASK)*.4;
imshow(M,[]);
%% view registration between movies
for e = 1:numel(regImage)
    imshow(regImage{e},[]);
    title(num2str(e))
    drawnow
    waitforbuttonpress
end
%% get mask for loading 
for e = 1:numel(regImage)
    miniStack(:,:,:,e) = regImage{e};
end
KK = mean(miniStack,4);
imshow(KK,[]);
%imshow(MASK,[]);
h = impoly();
pos = h.getPosition;
%%
M = poly2mask(pos(:,1),pos(:,2),size(KK,1),size(KK,2));
imshow(M);
sM = imresize(M,.3);
%% view sim and smoothed raw together
close all
samfret = 1;
linePlot = figure;
moviePlot = figure;
for sam = 1:size(SIM,5)
    SKIP = 20;
    for tm = 1:SKIP:size(SIM,4)
        figure(moviePlot);
        simD = SIM(:,:,:,tm,sam);
        rawD = MOV(:,:,:,tm,sam);
        img = cat(1,simD,rawD);
        img = cat(3,img,zeros(size(img,1),size(img,2)));
        imshow(img,[])
        drawnow
    end
    if samfret
        figure(linePlot)
        pointSIM = squeeze(SIM(x,y,:,:,sam));
        pointRAW = squeeze(MOV(x,y,:,:,sam));
        pointSIM = pointSIM(2,:).*pointSIM(1,:).^-1;
        pointRAW = pointRAW(2,:).*pointRAW(1,:).^-1;
        plot(pointSIM,'r');
        hold on
        plot(pointRAW,'k');
        hold on
    end
end
%% load genotype data
G = xlsread('/home/nate/Downloads/flg22 RT information for movies.xlsx');
%%  hand tag genotypes
WT10 = [1 2 3 9 12 11 12 19 23];
%% look at PC values over time
close all
CL = {'r' ,'b','g'};
for p = 1:size(ttC,1)
    tmp = squeeze(ttC(p,:,:));
    figure;
    plot(tmp,CL{p});
    hold on
    plot(tmp(:,end),'k')
    tmp = mean(tmp,2);
    plot(tmp,CL{p},'LineWidth',3)
    hold all
end
%% sweep eig-movie
close all
utC = mean(ttC,3);
%utC(1,:) = 0;
%utC(3,:) = 0;
sweep = PCA_BKPROJ_T(utC,tE,tU);
sweep = reshape(sweep,[SZ(1:3) size(sweep,2)]);
for t = 1:size(sweep,4)
    fretS = sweep(:,:,2,t).*sweep(:,:,1,t).^-1;
    fretS = fretS.*sM;
    imshow(fretS,[1.2 2.4])
    title(num2str(t))
    drawnow
end
%%
%sMOV = MOV;
%sMOV = zeros([128 410 3 301 29]);
tmp = tmpM(:,:,:,1);
tmp = imresize(tmp,.3);
sMOV = zeros([size(tmp) size(MOV,4) size(MOV,5)]);
for s = 1:size(MOV,5)
    tmpM = MOV(:,:,:,:,s);
    G = [];
    parfor t = 1:size(MOV,4)
        tmp = tmpM(:,:,:,t);
        tmp = imresize(tmp,.3);
        tmp = bsxfun(@times,tmp,sM);
        G(:,:,:,t) = imfilter(tmp,fspecial('disk',5));
        %sMOV(:,:,:,t,s) = imfilter(sMOV(:,:,:,t,s),fspecial('disk',5));
        t
    end
    sMOV(:,:,:,:,s) = G;
    s
end
%%
U = mean(MOV,5);
%% view average movie
close all
FRET = 1;
for loop = 1:5
    for t = 1:5:size(U,4)
        img = U(:,:,:,t);
        img = cat(3,img,zeros(size(img,1),size(img,2)));
        img = bsxfun(@times,img,sM);
        %img = img .*M;
        img(img(:) > 2) = 0;
        img = imresize(img,.4^-1);
        if FRET 
            img = img(:,:,2).*img(:,:,1).^-1;
            img(isnan(img(:))) = 0;
            imshow(img,[.5 2]);
        else
            imshow(img)
        end
       
        title(num2str(t));
        drawnow
    end
end
%% get wild type mean
WTidx = [1 2 3 13:18];
wtU = mean(SIM(:,:,:,:,WTidx),5);
%%
RBODidx = [4:7 (7:12)+12];
rbodU = mean(SIM(:,:,:,:,RBODidx),5);
%%
ACAidx = [8:12 (12:17)+12];
acaU = mean(SIM(:,:,:,:,ACAidx),5);
%% make fet for both
mutFRET = squeeze(rbodU(:,:,2,:).*rbodU(:,:,1,:).^-1);
acaFRET = squeeze(acaU(:,:,2,:).*acaU(:,:,1,:).^-1);
wtFRET = squeeze(wtU(:,:,2,:).*wtU(:,:,1,:).^-1);
%%
sM = imresize(M,.3);
%% probe spot
close all
tmp = reshape(tU,SZ(1:3));
[y x v] = impixel(tmp(:,:,1));
figure
sam = 1;
pointSIM = squeeze(SIM(x,y,:,:,sam));
pointRAW = squeeze(MOV(x,y,:,:,sam));
pointSIM = pointSIM(2,:).*pointSIM(1,:).^-1;
pointRAW = pointRAW(2,:).*pointRAW(1,:).^-1;
plot(pointSIM,'r')
hold on
plot(pointRAW,'k')
%% view total fret
close all
for loop = 1:10
    for e = 1:2:size(wtFRET,3)
        fr1 = wtFRET(:,:,e);
        %fr1 = imfilter(fr1,fspecial('disk',21));
        fr1(fr1(:) > 2) = 0;
        fr1(isnan(fr1(:))) = 0;


        fr2 = acaFRET(:,:,e);
        %fr2 = imfilter(fr2,fspecial('disk',21));
        fr2(fr2(:) > 2) = 0;
        fr2(isnan(fr2(:))) = 0;

        fr3 = mutFRET(:,:,e);
        %fr3 = imfilter(fr3,fspecial('disk',21));
        fr3(fr3(:) > 2) = 0;
        fr3(isnan(fr3(:))) = 0;

        RGB = cat(3,fr1,fr2,fr3);
        RGB = bsxfun(@times,RGB,sM);
        RGB = RGB -1;
        RGB = imresize(RGB,.3^-1);
        imshow(RGB,[]);
        title(num2str(e))
        drawnow

    end
end
%% make eig mon
close all
VEC = [];
myU = reshape(tU,SZ(1:3));
myU = myU(:,:,2).*myU(:,:,1).^-1;
%VEC = myU.*sM;
VEC = myU;
VEC = bindVec(VEC);
mag = [1 1 1 1 1];
for p = 1:3
    vec = reshape(tE(:,p),SZ(1:3));
    vec = vec(:,:,2).*vec(:,:,1).^-1;
    %vec = vec.*sM;
    vec(find(sM)) = bindVec(vec(find(sM)));
    vec = vec*mag(p);
    %vec = cat(3,vec,zeros(size(vec,1),size(vec,2)));
    %vec = bindVec(vec);
    vec = vec.*sM;
    VEC = [VEC;vec];
end
imshow(VEC,[]);
%%
%%

%%

%%
close all
for tr = 1:size(SIM,5)
    for f = 1:size(SIM,4)
        img = SIM(:,:,:,f,tr);
        img = cat(3,img,zeros(size(img,1),size(img,2)));
        imshow(imresize(img,.3^-1),[]);
        title(num2str(f));
        drawnow
    end
end
%% plot movie strands
%close all
IDX{1} = [1 2 3 13:18];
IDX{2} = [4:7 (7:12)+12];
IDX{3} = [8:12 (12:17)+12];
CL = {'k' 'g' 'r'};
h1 = figure;
h2 = figure;
hN{1} = figure;
hN{2} = figure;
hN{3} = figure;
hN{4} = figure;
hN{5} = figure;
mag = [1 1 1 1 1];
for g = 1:numel(IDX)
    for tr = 1:numel(IDX{g})
        figure(h1);
        plot(tC(1,:,IDX{g}(tr)),tC(2,:,IDX{g}(tr)),CL{g})
        hold on
        plot(tC(1,1,IDX{g}(tr)),tC(2,1,IDX{g}(tr)),['b*'])
        plot(tC(1,28:32,IDX{g}(tr)),tC(2,28:32,IDX{g}(tr)),['c.'])
        hold on
        figure(h2)
        plot(tC(3,:,IDX{g}(tr)),CL{g});
        hold on
    end
    for p = 1:5
        figure(hN{p});
        hold on
        Utr{g} = mag(p)*squeeze(mean(tC(p,:,IDX{g}),3));
        Str{g} = squeeze(std(tC(p,:,IDX{g}),1,3)*numel(IDX{g})^-.5);
        errorbar(Utr{g},Str{g},CL{g});
    end
end
%%

%%
for e = 1:size(tC,3)
    plot(tC(1,:,e),tC(2,:,e),'r')
    hold on
end
%%
%%
close all
for t = 1:size(U,4)
    %img = U(:,:,2,t);%.*U(:,:,2,t).^-1;
    img = U1(:,:,t).^-1;%.*U(:,:,2,t).^-1;
    img = img .*M;
    img(img(:) > 2) = 0;
    imshow(img,[.3 1.2]);
    title(num2str(t));
    drawnow
end
%%
close all
img = U(:,:,1,:).*U(:,:,2,:).^-1;
plot(squeeze(img(156,704,:,:)))
%%
close all
img = U1;
plot(squeeze(img(156,704,:)))
%% look at each image reg to first movie
close all
for e = 1:numel(regImage)
    U1 = mean(regMovie{1},3);
    U2 = regImage{e};
    RGB = cat(3,U1,U2,zeros(size(U1)));
    imshow(RGB,[]);
    drawnow
    waitforbuttonpress
end
%%
imageFile = '/mnt/tetra/nate/12.6.16.mdb/120616spl7.lsm';
info = imfinfo(imageFile);
clear I
cnt = 1;
for t = 1:2:numel(info)
    tmp = imread(imageFile,[],t);
    imshow(tmp,[]);
    I(:,:,:,cnt) = tmp;
    cnt = cnt + 1;
    drawnow
end

%% 
for t = 1:size(I,4)
    imshow(I(:,:,1,t),[])
    drawnow
end
%% 
for t = 1:20
    B = I(:,:,1,t) > 20;
    imshow(B,[])
    drawnow
end
%%
I = double(I);
%%
for t = 1%:size(I,4)
    temp = im2colF(I(:,:,1,t),[41 41], [15 15]);
end
%%
auto = trainAutoencoder(temp,3);
%% nate test
temp = im2colF(I(:,:,1,t),[41 41], [1 1]);
dec = encode(auto,temp);
decI(:,:,1) = col2im(dec(1,:),[41 41],[size(I,1) size(I,2)]);
decI(:,:,2) = col2im(dec(2,:),[41 41],[size(I,1) size(I,2)]);
decI(:,:,3) = col2im(dec(3,:),[41 41],[size(I,1) size(I,2)]);
%% align images
[optimizer,metric] = imregconfig('monomodal');
optimizer = registration.optimizer.RegularStepGradientDescent();
optimizer.MaximumIterations = 500;
%metric = registration.metric.MattesMutualInformation();
TAGNUM = numel(info)-1;
moving_reg = [];
parfor e = 1:size(I,4)
    idx = (e-1)*2+1;
    %target = rgb2gray(double(imread(imageFile,[],TAGNUM))/255);
    %source = rgb2gray(double(imread(imageFile,[],idx))/255);
    target = double(imread(imageFile,[],TAGNUM))/255;
    source = double(imread(imageFile,[],idx))/255;
    source = source(:,:,1);
    target = target(:,:,1);
    
    %imshow(source,[]);
    %figure;imshow(target);
    moving_reg(:,:,e) = imregister(source,target,'translation',optimizer,metric,'DisplayOptimization',1);
    e
end
%% look at reg movie
close all
for e = 1:numel(moving_reg)
    target = double(imread(imageFile,[],TAGNUM))/255;
    
    toShow = cat(3,target(:,:,1),moving_reg(:,:,e),zeros(size(target,1),size(target,2)));
    imshow(toShow,[]);
    drawnow
end
%% look at non reg movie
close all
for e = 1:size(I,4)
    imshow(I(:,:,:,e)/255,[]);
    drawnow
end
%% 