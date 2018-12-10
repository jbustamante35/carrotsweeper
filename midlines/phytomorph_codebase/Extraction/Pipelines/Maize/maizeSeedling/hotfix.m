
FilePath = '/home/nate/Downloads/funnyMN/';
FileList = {};
FileExt = {'nef'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
STACK = [];
IS = [];
newSZ = [3280 4992];
for e = 1:50%numel(FileList)
    try
        I = imread(FileList{e});
        I = imresize(I,newSZ);
        I(1:1400,:,:) = [];
        IS(:,:,:,e) = I;
        I = imresize(I,.25);
        sz = size(I);
        I = reshape(I,[prod(sz(1:2)) sz(3)]);
        
        STACK = [STACK;I];
        
        e
        size(IS)
        rm(e) = false;
    catch
        rm(e) = true;
    end
end
%% dev JUNK
qrCropBox = {};
BUMP = 50;
uS2 = [];
sS2 = [];
uS1 = [];
sS1 = [];
toLook = 2;
for e = 1:numel(FileList)
    I = imread(FileList{e});
    
    
    
    
    I = imresize(I,newSZ);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [msg,qrCropBox{e}] = getQRcode(I);
    TOPloc = qrCropBox{e}(2) + qrCropBox{e}(4);
    TOPloc = TOPloc + BUMP;
    
    
    oI = I;
    
    I(1:TOPloc,:,:) = [];
    
    
    tmpHSV = rgb2hsv(I);
    
    %uS2(:,e) = mean(tmpHSV(:,:,toLook),2);
    %sS2(:,e) = std(tmpHSV(:,:,toLook),1,2);
    uS1(:,e) = mean(tmpHSV(:,:,toLook),1)';
   
    sS1(:,e) = std(tmpHSV(:,:,toLook),1,1)';
    e
    sig = sS1(:,e).*uS1(:,e).^-1;
    ssig = imfilter(sig,fspecial('average',[501 1]),'circular');
    sv = sort(ssig);
    
   
    ssig(isnan(ssig)) = min(ssig(:));
    BV = zeros(1,size(I,2));
   
    
    ssig = bindVec(ssig);
    
    
    G = rgb2gray(I);
    uG = mean(imcomplement(G),1);
    uG = imfilter(uG,fspecial('average',[501 1]),'circular');
    ssig = (ssig'.*uG)';
    ssig = bindVec(ssig);
    threshV = graythresh(ssig);
    BV = (ssig < threshV)';
    
    
    
    
    BV = logical([zeros(size(BV));BV;zeros(size(BV))]);
    BV = imclearborder(BV);
    R = regionprops(logical(BV),'Centroid');
    BV = BV(2,:);
    BV = zeros(size(BV));
    for e = 1:numel(R)
        BV(1,round(R(e).Centroid(1))) = 1;
    end
    BV(1) = 1;
    BV(end) = 1;
    
    
    pts = find(BV);
    
    for b = 1:(numel(pts)-1)
        CB{b} = [pts(b) TOPloc pts(b+1)-pts(b) size(oI,2)-TOPloc];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for b = 1:numel(CB)
        J = imcrop(oI,CB{b});
        
        
        sJ = imfilter(J,fspecial('disk',50),'replicate');
        tmpHSV = rgb2hsv(sJ);
        tmpU = mean(tmpHSV(:,:,2),2);
        tmpS = std(tmpHSV(:,:,2),1,2);
        tmpSIG = tmpU.*tmpS;
        tmpSIG = bindVec(tmpSIG);
        tmpBV = tmpSIG > graythresh(tmpSIG);
        SUB = sum(tmpBV);
        CB{b}(4) = CB{b}(4) - SUB;
        
        
    end
    
    
    
    
    
    
    BV = logical(repmat(BV,[size(oI,1) 1]));
    out = flattenMaskOverlay(oI,BV);
    imshow(out,[]);
    hold on
    plot(1:size(I,2),ones(1,size(I,2))*TOPloc,'r')
    plot(1:size(I,2),ssig*1000,'g')
    plot(1:size(I,2),uS1(:,e)'*10000,'m')
    hold off
    drawnow
    
end
%%
close all
for e = 1:size(uS2,2)
    sig = uS2(:,e).*sS2(:,e);
    I = imread(FileList{1});
    [msg,qrCropBox] = getQRcode(I);
    
    baseLine = imfilter(sS2(:,e),fspecial('average',[201 1]),'replicate');
    nor = baseLine-sS2(:,e);
    TH = abs(nor);
    plot(sS2(:,e),'k');
    hold on
    plot(baseLine,'r')
    
    plot(TH,'b')
    plot(sig,'g');
    hold off
    %waitforbuttonpress
    pause(.2)
end
%% 

I = norSTACK(:,:,:,e);

[mask] = cropImages_ver2(I,950,TOPTHRESH,GMModels{3},level2,levelA);


%%


[returnI,boundingBoxes,connectedMASK,MASK,SKELETON] = cropImages_ver2(FileList{1},950,GMModels{3},level2,levelA);
%%
CalliList{1} = '/home/nate/Downloads/{Plot_1013}{Experiment_B97RILS}{Planted_1-3-18}{GenoA_Z001E0126}{GenoB_Z001E0163}{GenoC_Z001E0033}{Treatment_Control}{PictureDay_4}.nef';
CalliList{1} = '/home/nate/Downloads/{Plot_1013}{Experiment_B97RILS}{Planted_1-3-18}{GenoA_Z001E0126}{GenoB_Z001E0163}{GenoC_Z001E0033}{Treatment_Control}{PictureDay_4}.nef';
%%
otherList{1} = '/home/nate/Downloads/{Plot_2006}{Experiment_73}{Planted_12-8-2017}{SeedSource_DI1999-1}{SeedYear_2015}{Genotype_B73}{Treatment_Control}{PictureDay_15}.nef';
otherList{1} = '/home/nate/Downloads/seedlingTest/{Plot_2007}{Experiment_73}{Planted_12-8-2017}{SeedSource_DI1999-1}{SeedYear_2015}{Genotype_B73}{Treatment_Control}{PictureDay_15}.tiff';
otherList{1} = '/home/nate/Downloads/seedlingTest/{Plot_1307}{Experiment_57}{Planted_8-8-2017}{SeedSource_DI2125-4}{SeedYear_2015}{Genotype_Ki11}{Treatment_W_1_2d_cold}{PictureDay_9}.tiff';

smartMain_v4(otherList{1},'./output/',[],GMModels{3},level2,levelA,levelD,kE,kU);
%%
for e = 1:numel(FileList)
    smartMain_v4(FileList{e},'./output/',[],GMModels{3},level2,levelA,levelD,kE,kU);
end
%% last cut
KK = [];
sigSZ = 21;
KK2 = zeros(2000000,(sigSZ*2+1)*3);

cnt = 1;
for e = 1:numel(FileList)
    [I,returnI,boundingBoxes,connectedMASK,MASK,SKELETON] = cropImages_ver3(FileList{1},950,GMModels{3},level2,levelA,[]);
    %{
    for k = 1:size(returnI)
        [c r v] = impixel(returnI{e});
        stripV = returnI{e}(r,c-sigSZ:c+sigSZ,:);
    end
%}
    
    %{
    for e = 1:3
        
        stmpI = imfilter(returnI{e},fspecial('average',[1 30]),'replicate');
        [newMask] = applyClusterStage(stmpI,MASK{e},levelD);
    end
    %}

    for o = 1:3
        %{
        HSVGG = rgb2hsv(returnI{o});
        [mask] = gatherPlantnHoods(returnI{o},MASK{o},[0 21],kE,kU,levelD,[1:3],levelDD);
        %}
        
        

        kI = [];
        idx = find(MASK{o});
        [ri ci] = find(MASK{o});
        toSample = padarray(returnI{o},[sigSZ sigSZ],'replicate','both');
        for p = 1:numel(ri)
            if cnt < size(KK2,1)
                stVec = toSample(ri(p)+sigSZ,(ci(p)+sigSZ-sigSZ):(ci(p)+sigSZ+sigSZ),:);
                stVec = sort(stVec,2);

                KK2(cnt,:) = stVec(:);

                cnt = cnt + 1;
            end
                %KK2 = cat(1,KK2,stVec(:)');
        end
        
        o
        cnt
        
        %{
        stmpI = imfilter(returnI{o},fspecial('average',[1 30]),'replicate');
        for k = 1:3
            tmp = stmpI(:,:,k);

            

            kI = [kI tmp(idx)];
        end
        KK = [KK;kI];
%}
    end

end
%%
FFF = zeros(size(KK2,1),9);
parfor e = 1:size(KK2,1)
    tmp = KK2(e,:);
    tmp = reshape(tmp,[1 2*sigSZ+1 3]);
    f1 = mean(tmp(1,(end-4):end,:)) - mean(tmp(1,1:5,:));
    f2 = mean(tmp,2);
    f4 = std(tmp,1,2);
    f3 = [f1(:) f2(:) f4(:)];
    FFF(e,:) = f3(:)';
e
end
%%
zscore(FF);
%%
[kS kC kU kE kL kERR kLAM] = PCA_FIT_FULL(KK2,3);
%%
K = 3;
SAM = 10;
options = statset('Display','iter','MaxIter',500);
%levelD = fitgmdist(kC(:,[1]),K,'Start','plus','Options',options);
levelD = fitgmdist(kC(:,1:3),K,'Start','plus','Options',options,'RegularizationValue',.000001);
%%
% sub cluster BR
K = 4;
tempyIDX = cluster(levelD,FFF(:,1:3));
SAM = 10;
options = statset('Display','iter');
levelDD = fitgmdist(FFF(tempyIDX==3,4:6),K,'Start','plus','Options',options);

%%
close all
for e = 1:size(kC,2)
    ksdensity(kC(1:5:end,e))
    drawnow

end

%%
IS = IS(:,:,:,1:2:end);
%% correct brightness
uS = squeeze(mean(IS,4)/255);
uHSV = rgb2hsv(uS);
H = imhist(uHSV(:,:,3));
CS = [];
CS2 = [];
CS3 = [];
for k = 1:3
    KH{k} = imhist(uS(:,:,k));
end



parfor e = 1:size(IS,4)
    
    
    toOp = IS(:,:,:,e)/255;
    toOp(:,end-50:end,:) = [];
    
    RGB = toOp;
    HSV = rgb2hsv(RGB);
    
    fRGB = imfilter(RGB,fspecial('average',[1 31]),'replicate');
    %{
    smoothed = imresize(toOp,.25);
    smoothed = imfilter(smoothed,fspecial('disk',25),'replicate');
    smoothed = imresize(smoothed,4);
    
    
    
    %STORE = [];
    inI = toOp;
    for k = 1:1
        HSV = rgb2hsv(inI);
        %J = imhistmatch(HSV(:,:,3),uHSV(:,:,3),255);
        %HSV(:,:,3) = J;
        HSV(:,:,3) = mean(mean(uHSV(:,:,3)))*ones(size(HSV,1),size(HSV,2));
        RGB = hsv2rgb(HSV);
        %RGB = mean(cat(4,RGB,toOp),4);
        inI = RGB;
        %STORE(:,:,:,k) = inI;
        k
    end
    
    RGB = imfilter(RGB,fspecial('disk',11),'replicate');
    %}
    
    
    
    %{
    RGB = imhistmatch(RGB,uS,255);
    %{
    for k = 1:3
        %RGB(:,:,k) = histeq(RGB(:,:,k),KH{k});
        RGB(:,:,k) = imhistmatchn(RGB(:,:,k),uS(:,:,k));
    end
    %}
    
    HSV = rgb2hsv(RGB);
    %{
    DIS(:,:,1) = sum(RGB.^2,3).^.5;
    gRGB = RGB;
    gRGB(:,:,2) = 1 - RGB(:,:,2);
    DIS(:,:,2) = sum(RGB.^2,3).^.5;
    gRGB = RGB;
    gRGB(:,:,3) = 1 - RGB(:,:,3);
    DIS(:,:,3) = sum(RGB.^2,3).^.5;
    %}
    
    out = cat(2,RGB,IS(:,:,:,e)/255);
    
    %imshow(out,[]);
    %}
    CS2(:,:,:,e) = HSV;
    CS(:,:,:,e) = RGB;
    CS4(:,:,:,e) = fRGB;
    %CS3(:,:,:,e) = smoothed;
    %drawnow
    e
end
%%
smHSV = permute(CS2,[1 2 4 3]);
sz = size(smHSV);
smHSV = reshape(smHSV,[prod(sz(1:3)) sz(4)]);

%%
smF = permute(CS4,[1 2 4 3]);
sz = size(smF);
smF = reshape(smF,[prod(sz(1:3)) sz(4)]);
%%
smS = permute(CS3,[1 2 4 3]);
sz = size(smS);
smS = reshape(smS,[prod(sz(1:3)) sz(4)]);

%%
smRGB = permute(CS,[1 2 4 3]);
sz = size(smRGB);
smRGB = reshape(smRGB,[prod(sz(1:3)) sz(4)]);
%%
SAM = 100;
for k = 1:4
    GMModels{k} = fitgmdist([smRGB(1:SAM:end,:),smHSV(1:SAM:end,[2 3])],k,'Start','plus','Options',options);
    AIC(k)= GMModels{k}.AIC;
end
%%

%%
%sleM1 = GMModels{3};
%idxT = cluster(sleM1,[smRGB,smHSV(:,[2 3])]);
fidx1 = find(idxT==3);
numel(fidx1)
SAM = 10;
K = 4;
options = statset('Display','iter');
%level2 = fitgmdist([smRGB(fidx1(1:SAM:end),:),smHSV(fidx1(1:SAM:end),[2 3])],K,'Start','plus','Options',options,'RegularizationValue',.01);
level2 = fitgmdist([smF(fidx1(1:SAM:end),:)],K,'Start','plus','Options',options,'RegularizationValue',.000001);
%%
K=3;
SAM = 100;
options = statset('Display','iter');
levelA = fitgmdist(A,K,'Start','plus','Options',options);


%%
fidx1 = find(idxT==3);
fidx3 = find(idxT2==1);
numel(fidx3)
K = 3;
SAM = 100;
options = statset('Display','iter');
level4 = fitgmdist([smF(fidx1(fidx3),:)],K,'Start','plus','Options',options,'RegularizationValue',.00001);


%%
fidx1 = find(idxT==3);
idxT2 = cluster(level2,[smRGB(fidx1,:),smHSV(fidx1,[2 3])]);
fidx2 = find(idxT2==2);
numel(fidx2)
K = 3;
SAM = 100;
options = statset('Display','iter');
level3 = fitgmdist([smF(fidx1(fidx2(1:SAM:end)),:)],K,'Start','plus','Options',options,'RegularizationValue',.00001);




%%
SAM = 100;
K = 3;
MAJOR = fitgmdist(smS(1:SAM:end,:),K,'Start','plus','Options',options);
%IDX_MAJOR = cluster(MAJOR,smS);
%%
SAM = 100;
K = 4;
minusBright = fitgmdist(smRGB(1:SAM:end,:),K,'Start','plus','Options',options);
%IDX_MAJOR = cluster(MAJOR,smS);
%%

%%
SAM = 100;
K = 2;
BR = fitgmdist(CS2(1:SAM:end,3),K,'Start','plus','Options',options);
IDX_BR = cluster(BR,CS2(:,3));


%%
% sub cluster BR
K = 3;
sum(IDX_BR==1)
fidx1 = find(IDX_BR==1);
SAM = 10;
options = statset('Display','iter');
CL_sub1 = fitgmdist(CS(fidx1(1:SAM:end),:),K,'Start','plus','Options',options);


% sub cluster BR
K = 3;
sum(IDX_BR==2)
fidx1 = find(IDX_BR==2);
SAM = 10;
options = statset('Display','iter');
CL_sub2 = fitgmdist(CS(fidx1(1:SAM:end),:),K,'Start','plus','Options',options);

%%
close all
norSTACK = CS;

uS = mean(norSTACK,4);
norHSV = rgb2hsv(uS);
norBright = norHSV(:,:,3);
%BRF = max(CS2(:,:,3),[],4);
%A = [];
for e = 1:50
    toMask = norSTACK(:,:,:,e);
    [mask] = generatePlantMasks_ver2(toMask,GMModels{3},level2,levelA);
    out = flattenMaskOverlay(toMask,mask);
    imshow(out,[]);
    drawnow
    
    %{
    
    
    toFix = norSTACK(:,:,:,e);
    fRGB = imfilter(toFix,fspecial('average',[1 31]),'replicate');
    toFixHSV = rgb2hsv(toFix);
    
    
    sz = size(toFix);
    Q1 = reshape(toFix,[prod(sz(1:2)) sz(3)]);
    sz = size(toFixHSV);
    Q2 = reshape(toFixHSV,[prod(sz(1:2)) sz(3)]);
    sz = size(fRGB);
    Q3 = reshape(fRGB,[prod(sz(1:2)) sz(3)]);
    [M_idx,nlogl,P] = cluster(GMModels{3},[Q1 Q2(:,[2 3])]);
    M_idx = reshape(M_idx,sz(1:2));
    fidx = find(M_idx==3);
    
    %[M2_idx] = cluster(level2,[Q1(fidx,:) Q2(fidx,[2 3])]);
    [M2_idx] = cluster(level2,[Q3(fidx,:)]);
    
    %{
    gidx = find(M2_idx==2);
    [M3_idx] = cluster(level3,Q3(fidx(gidx),:));
    
    gidx2 = find(M2_idx==1);
    [M3_idx2] = cluster(level4,Q3(fidx(gidx2),:));
    %}
    MT = M_idx;
    MT(fidx) = M2_idx+3;
    mask = MT == 7 | MT == 5;
    
    
    
    qidx = find(mask);
    [M_idxA] = cluster(levelA,[Q1(qidx,:)]);
    mask = double(mask);
    mask(qidx) = M_idxA;
    base = mask == 1;
    
    base = bwareaopen(base,50);
    bidx = find(base);
    
    ridx = regionprops(imdilate(mask==2,strel('disk',1)),'PixelIdxList');
    
    for e = 1:numel(ridx)
        if ~isempty(intersect(ridx(e).PixelIdxList,bidx))
            base(ridx(e).PixelIdxList) = 1;
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    mask = bwareaopen(mask,100);
    %mask = imclose(mask,strel('disk',11,0));
    %mask = imclearborder(mask);
    
    %mask = bwlarge(mask,3);
    
    
    %{
    MT(fidx(gidx)) = M3_idx+6;
    %MT(fidx(gidx2)) = M3_idx2+9;
    %}
    %mask = MT == 7 | MT == 8;
    out = flattenMaskOverlay(toFix,base);
    %{
    lidx = find(mask);
    A = [A;Q1(lidx,:)];
    %}
    
    imshow(out,[]);
    drawnow
    
    %}
    %{
    MTL = label2rgb(MT);
    
    smoothed = imresize(toFix,.25);
    smoothed = imfilter(smoothed,fspecial('disk',5),'replicate');
    smoothed = imresize(smoothed,4);
    
    HSV = rgb2hsv(toFix);
    tmpV = HSV(:,:,3);
    fidx = find(tmpV > .5);
    tmpV(fidx) = .7;
    HSV(:,:,3) = tmpV;
    
    
    RGB = hsv2rgb(HSV);
    
    sz = size(RGB);
    RGB = reshape(RGB,[prod(sz(1:2)) sz(3)]);
    [rgb_idx,nlogl,rgb_P] = cluster(minusBright,RGB);
    rgb_P = reshape(rgb_P,[sz(1:2) size(rgb_P,2)]);
    rgb_idx = reshape(rgb_idx,sz(1:2));
    rgb_label = label2rgb(rgb_idx);
    
    
    sz = size(smoothed);
    smoothed = reshape(smoothed,[prod(sz(1:2)) sz(3)]);
    %M_idx = cluster(MAJOR,smoothed);
    [M_idx,nlogl,P] = cluster(MAJOR,smoothed);
    
    P = reshape(P,sz);
    M_idx = reshape(M_idx,sz(1:2));
    
    
    [fixed D] = fixImage(toFix,norBright,norSTACK(:,:,:,1),10);
    norSTACK(:,:,:,e) = fixed;
    %}
end

%% try cov method
func = @(X)colorCOV(X,1);
W = colfilt(CS(:,:,:,1),[11 11 3],'sliding',func);
%%
close all
for e = 1:size(DIFF,3)
    imshow(DIFF(:,:,e),[min(DIFF(:)) max(DIFF(:))]);
    drawnow
    pause(1)
end

%%




HSV = rgb2hsv(RGB);
tmp = HSV(:,:,3);
idx1 = cluster(BR,tmp(:));

fidx1 = find(idx1==1);
fidx2 = find(idx1==2);

sz = size(RGB);
rgb = reshape(RGB,[prod(sz(1:2)) sz(3)]);
idx2 = cluster(CL_sub1,rgb(fidx1,:));
idx3 = cluster(CL_sub2,rgb(fidx2,:));


IDXM = zeros(size(idx1));

IDXM(fidx1) = idx2;
IDXM(fidx2) = 3+idx3;
IDXM = reshape(IDXM,[sz(1:2)]);
RGBL = label2rgb(IDXM);


%% reshape CS
CS = permute(CS,[1 2 4 3]);
sz = size(CS);
CS = reshape(CS,[prod(sz(1:3)) sz(4)]);
%%
STACK = CS;
%%
STACK = double(STACK);

%%
SAM = 10;
K = 3;
options = statset('Display','iter');
GMModel = fitgmdist(STACK(1:SAM:end,:),K,'Start','plus','Options',options);
%%
K = 3;
IDX = cluster(GMModel,STACK);
kidx = (IDX == 1);
options = statset('Display','iter');
GMModel2 = fitgmdist(STACK(kidx,:),K,'Start','plus','Options',options);
%%
close all
for e = 1:numel(FileList)
    I = imread(FileList{e});
    
    I = imresize(I,newSZ);
    I(1:1400,:,:) = [];
    
    HSV = rgb2hsv(double(I)/255);
    J = histeq(HSV(:,:,3),H);
    HSV(:,:,3) = J;
    RGB = hsv2rgb(HSV);
    for k = 1:3
        RGB(:,:,k) = histeq(RGB(:,:,k),KH{k});
    end
    I = RGB;
    
    sz = size(I);
    tmpI = reshape(I,[prod(sz(1:2)) sz(3)]);
    idx = cluster(GMModel,double(tmpI));
    idx = reshape(idx,sz(1:2));
    RGB = label2rgb(idx);
    mask = idx == 1;
    mask = bwareaopen(mask,1000);
    %mask = imclearborder(mask);
    
    midx = find(mask);
    idx2 = cluster(GMModel2,double(tmpI(midx,:)));
    mask2 = zeros(size(mask));
    mask2(midx) = idx2 == 3;
    
    out = flattenMaskOverlay(double(I)/255,mask);
    out = flattenMaskOverlay(out,logical(mask2),.5,'b');
    imshow(out,[]);
    drawnow
end
%%
for e = 1:size(IS,4)
    I = IS(:,:,:,e);
    [r{e} c{e} v{e}] = impixel(I/255);
    
end
%%
