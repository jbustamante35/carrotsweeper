%% 0) look for images
% arabidopsis
FilePath =  '/mnt/piranhapuffer1/sscam/Root Measure/High Res/Baseline Images/';
% tomato
%FilePath = '/mnt/piranhapuffer2/Bessie/Tomato ILs/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% 1) sort SET
for s = 1:numel(SET)
    NM = [];
    for e = 1:numel(SET{s})
        [pth nm ext] = fileparts(SET{s}{e});
        NM(e) = str2num(nm);        
    end
    [~,sidx] = sort(NM);
    SET{s} = SET{s}(sidx);
end
%% JUNK TESTT
img = imread(SET{6}{1});
extractEdgeCurves(img);
eC(SET{6}{1})
%% learning manifold
lm = lManifold();
lm.setinitComp(5);
lm.setGroupN(9);
lm.setmodelCompX(5);
lm.setmodelCompY(3);
T = goT();
T.setNhoodRho(30);
T.setNhoodRad(pi);
T.setNhoodDensity([30 100]);
T.generateH();
samRho = 20;
for e = 1:5%numel(SET)
    % read in the first image
    I = imread(SET{e}{1});
    % gather path
    imshow(I,[]);
    [x y V] = impixel();    
    close all
    % create path and tracer object
    curve = [x y]';
    curve = goT.reparameterize(curve);
    TB = goT.generateTangentBundle(curve);
    T.setImage(double(I));
    T.position = curve;
    T.direction = TB;
    % show the image and the frame bundle
    imshow(I,[])
    hold on;
    T.plotFrameBundle();
    samCurve = T.sampleCurveAtCurve(samRho);
    samCurve = reshape(samCurve,[size(samCurve,1)*size(samCurve,2) size(samCurve,3)]);
    samImage = T.sampleImageAtCurve();
    % take minus rho
    samImage = samImage(:,1:end-samRho);
    lm.addXY(samImage',samCurve');
end
%%
lm.cleanData();
lm.learn();
Yi = lm.predict(lm.rawX(10:100,:));


%% 2) load all first images
I = [];
S = [];
figure;
for i = 1:numel(SET)
    I = imread(SET{i}{1});
    
    if size(I,1) == 1040
        S = cat(3,S,I);
        imshow(I,[])
        drawnow
    end
end
%% 3) get model
SIG = zeros(size(S,1),1);
SIG(1) = 1;
for e  = 1:50%size(S,3)
    % get image
    img = S(:,:,e);
    img = double(img)/255;
    img(1,:) = [];
    % handle flip
    img = handleFLIP(img,[]);
    % extract along curve
    [curve{e} sample(:,:,:,e)] = extractCurves(img);
    e
end
% create the model for the tip
model = mean(sample,4);
%% view model
close all
for e = 1:size(model,1)
    imshow(squeeze(model(e,:,:)),[]);
    drawnow
end
%% extract curve and labels via model
SIG = zeros(size(S,1),1);
SIG(1) = 1;
SIG  = imdilate(SIG,ones(5,1));
PAD = 50;
c = 1;
for e  = 1:100%size(S,3)
    % get image
    img = S(:,:,e);
    img = double(img)/255;
    img(1,:) = [];
    % handle flip
    img = handleFLIP(img,[]);
    % extract along curve
    [curve{e},~,P] = extractCurves(img,model);
    % displace the length of the model
    P{1}{2} = P{1}{2} + (size(model,1)-1)/2;
    % get curve
    contour = fliplr(curve{e}(:));
    % extract midlines
    [db] = extractMidlines(contour,P);
end
%% 4) call stack process
oPath = '/mnt/spaldingdata/nate/devRuns/arabidopsis_run_2013.07.22/';
parfor s = 1:numel(SET)
    et = clock;
    stackP(SET{s},model,oPath);
    delta = clock - et;
end
%% other for rollaDex
SIG = zeros(size(S,1),1);
SIG(1) = 1;
SIG  = imdilate(SIG,ones(5,1));
PAD = 50;
c = 1;
for e  = 1:100%size(S,3)
    % get image
    img = S(:,:,e);
    img = double(img)/255;
    img(1,:) = [];
    % handle flip
    img = handleFLIP(img,[]);
    % extract along curve
    [curve{e} sample(:,:,:,e)] = extractCurves(img);
    
    %{
    
    % roll for X root on screen dim    
    level = graythresh(img);
    imgB = img < level;
    imgB = bwareaopen(imgB,2000);
    fidx = find(any(imgB,1));
    img = circshift(img,[0 -fidx(end)+250]);
    
    % align base
    tmp = mean(img(:,1:PAD),2);
    sig = cconv(-tmp,SIG,numel(tmp));
    sig = imfilter(sig,fspecial('average',21),'replicate');
    [~,sidx] = max(sig);
    
    img = circshift(img,[-sidx+300 0]);
    

    try
        H(:,:,c) = img;
        c = c + 1;
    end
    e
    % 
    
    
    imshow(img,[])
    drawnow
    %}
end
%% look at mean H
uH = mean(H,3);
imshow(uH,[])
%%
SIG = zeros(size(S,1),1);
SIG(1) = 1;
SIG  = imdilate(SIG,ones(5,1));
PAD = 50;
c = 1;
for e  = 1:size(S,3)
    % get image
    img = S(:,:,e);
    img = double(img)/255;
    img(1,:) = [];
    % handle flip
    img = handleFLIP(img,[]);
    % extract along curve
    [curve{e} sample{e}] = extractCurves(img);
end