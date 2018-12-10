%%% RUN STEP 1: CONFIRMED
%massDownload('/iplant/home/petersoapes/returns/return5/', '.mat','/mnt/tetra/nate/learnDNA/');
%massDownload('/iplant/home/petersoapes/returns/return5/', '.csv','/mnt/tetra/nate/learnDNA/');
%% dig for files
% run before publish code on Sept 21 2018
FilePath = '/mnt/tetra/nate/learnDNA/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% load loop
% needed from here - resize method to reshape the data
close all
newSZ = [500 61];
%run before publish code on Sept 21 2018
toLoad = numel(FileList);
%toLoad = 10000;
S = zeros([newSZ 3 toLoad]);
disp = false;
tm = [];
for e = 1:toLoad
    try
        tic
        
        load(FileList{e})
        % normalize vec
        for k = 1:size(vec,3)
            vec(:,:,k) = bindVec(vec(:,:,k));
        end
        vec = imresize(vec,newSZ);
        
        S(:,:,:,e) = vec;
        if disp
            imshow(vec,[]);
            drawnow
        end
        e
        rm(e) = false;
        tm(e) = toc;
        mean(tm)/60*(toLoad-e)
    catch
        rm(e) = true;
    end
end
%%%%%%%%%%%%%%%%%%
% view code
%% look at results from parloop load
%%%%%%%%%%%%%%%%%%
close all
for e = 1:size(S,4)
    imshow(S(:,:,:,e),[]);
    drawnow
    
end
%%
SBK = S;
%% remove bad images
%run before publish code on Sept 21 2018
rm = find(squeeze(any(any(any(isnan(S),1),2),3)));
S(:,:,:,rm) = [];
tmpList = FileList(1:toLoad);
tmpList(rm) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first level PCA
%run before publish code on Sept 21 2018
S = permute(S,[1 4 2 3]);
osz = size(S);
S = reshape(S,[prod(osz(1:2)) osz(3:4)]);
%% loop over colors
%run before publish code on Sept 21 2018
clear simS
clear C
clear U
clear E
numComponentsLevel1 = 3;
for e = 1:size(S,3)
    [simS(:,:,e),C(:,:,e),U{e},E{e}] = PCA_FIT_FULL(S(:,:,e),numComponentsLevel1);
end
%% reshape
simS = reshape(simS,osz);
simS = ipermute(simS,[1 4 2 3]);

csz = osz;
csz(3) = numComponentsLevel1;
simC = reshape(C,csz);
simC = ipermute(simC,[1 4 2 3]);

viewS = reshape(S,osz);
viewS = ipermute(viewS,[1 4 2 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%
% view code
% view the results of the first level PCA - stack next to the results raw
close all
for e = 1:size(viewS,4)
    tmp = [simS(:,:,:,e) viewS(:,:,:,e)];
    imshow(tmp,[]);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% second level PCA
simSr = permute(simS,[2 4 1 3]);
osz2 = size(simSr);
simSr = reshape(simSr,[prod(osz2(1:2)) osz2(3:4)]);
%% loop over colors level 2
clear simS2
clear C2
clear U2
clear E2 
numComponentsLevel2 = 3;
for e = 1:size(simS,3)
    [simS2(:,:,e),C2(:,:,e),U2{e},E2{e}] = PCA_FIT_FULL(simSr(:,:,e),numComponentsLevel2);
end
%% try to cluster data on second level PCA
noiseR = S;
szN = size(noiseR);
for e = 1:szN(4)
    mini(:,:,:,e) = imresize(noiseR(:,:,:,e),.20);
    e
end
miniSZ = size(mini);
mini = reshape(mini,[prod(miniSZ(1:3)) miniSZ(4)]);
[simMINI,CMINI,UMINI,MINI] = PCA_FIT_FULL_T(mini,3);
options = statset('Display','iter');
gm = fitgmdist(CMINI',3,'Options',options);
kidx = gm.cluster(CMINI');
%% reshape
simS2 = reshape(simS2,osz2);
simS2 = ipermute(simS2,[2 4 1 3]);

csz = osz2;
csz(3) = numComponentsLevel2;
simC2 = reshape(C2,csz);
simC2 = ipermute(simC2,[2 4 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%
% view code
% view the results of the first level PCA - stack next to the results raw
close all
for e = 1:size(viewS,4)
    tmp = [viewS(:,:,:,e) simS(:,:,:,e) simS2(:,:,:,e)];
    imshow(tmp,[]);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASK level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
vec = simC2(1,:,1,:);
ksdensity(vec(:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix S next -->? yes
S = reshape(S,osz);
S = ipermute(S,[1 4 2 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix S next -->? yes - run through twice
MS = S;
for e = 1:size(S,4)
    
    top = simC2(1,:,1,e);
    topMask = top > 0;
    
    
    topMask = repmat(topMask,[size(simS,1) 1]);

    MASKSTACK(:,:,e) = topMask;

    % need to define this before applying
    topMask = newM;

    MS(:,:,:,e) = bsxfun(@times,S(:,:,:,e),topMask);
    e
end
%%
newM = mean(MASKSTACK,3)>.8;
%% first level PCA MASK fix S next -->? yes
MS = permute(MS,[1 4 2 3]);
oszM = size(MS);
MS = reshape(MS,[prod(oszM(1:2)) oszM(3:4)]);
%% loop over colors fix S next -->? yes
clear simMS
clear MC
clear MU
clear ME 
numComponentsLevelMASK = 2;
for e = 1:size(MS,3)
    [simMS(:,:,e),MC(:,:,e),MU{e},ME{e}] = PCA_FIT_FULL(MS(:,:,e),numComponentsLevelMASK);
end
%% reshape MASK --> ? yes
simMS = reshape(simMS,oszM);
simMS = ipermute(simMS,[1 4 2 3]);

MS = reshape(MS,oszM);
MS = ipermute(MS,[1 4 2 3]);

csz = osz;
csz(3) = numComponentsLevelMASK;
simMC = reshape(MC,csz);
simMC = ipermute(simMC,[1 4 2 3]);
simMC = squeeze(simMC);
%%
%%%%%%%%%%%%%%%%%%
% view code
% view the results of the first level PCA - stack next to the results raw
close all
mag = 30
probe = 2;
fidx = find(kidx==probe);
%for i = 1:size(viewS,4)


initLabels = simMC(:,:,:,fidx);
initSZ = size(initLabels);
initLabels = reshape(initLabels,[initSZ(1) prod(initSZ(2:3)) initSZ(4)]);
initLabels = permute(initLabels,[2 1 3]);
initSZ = size(initLabels);
initLabels = reshape(initLabels,[prod(initSZ(1)) prod(initSZ(2:3))]);

options = statset('Display','iter');
gmL = fitgmdist(initLabels(:,1:10:end)',5,'Options',options);
lidx = gmL.cluster(initLabels');
lidx = reshape(lidx,initSZ(2:3));
%%
for i = 1:numel(fidx)
    
    tLAB = zeros(size(viewS,1),1);
    
    e = fidx(i);
    tmp = [viewS(:,:,:,e) simS(:,:,:,e) simS2(:,:,:,e) simMS(:,:,:,e) MS(:,:,:,e)];


    imshow(tmp,[]);
    hold on

    sortedRed = sort(simMC(:,1,1,e));
    threshold = mean(sortedRed(1:50));
    
    
    
    sig = bindVec(simMC(:,1,1,e));
    threshold = graythresh(sig);
    CH = sig > threshold;
    CH = bwlarge(CH);

    CHred = simMC(:,1,1,e);
    CHred = CHred(find(CH));
    

    CHgreen = simMC(:,1,2,e);
    gfidx = find(CH);
    CHgreen = CHgreen(find(CH));
    CHgreen = bindVec(CHgreen);
    CHgreen = CHgreen - mean(CHgreen);
    gSIG = bindVec(CHgreen);
    
    gidx = find(gSIG > graythresh(gSIG));
    
    %CHgreen = CHgreen*mean(CHred);

    dG = zeros(size(sortedRed));
    dG(find(CH)) = CHgreen;





    plot(simMC(:,1,1,e)*mag + 61*2,1:500,'r')
    plot(simMC(:,2,1,e)*mag*2 + 61*2,1:500,'m')
    
    plot(simMC(:,1,2,e)*mag + 61*2,1:500,'g')
    plot(simMC(:,2,2,e)*mag*2 + 61*2,1:500,'y')
    
    plot(simMC(:,1,3,e)*mag + 61*2,1:500,'b')
    plot(simMC(:,2,3,e)*mag + 61*2,1:500,'c')
    plot(CH*mag + 20,1:500,'r')
    plot(-dG*2*mag + 20,1:500,'g')

    
    
    [~,CENTRO(i)] = max(simMC(:,1,3,e));
    
    
    
    tLAB(find(CH)) = 1;
    tLAB(CENTRO(i)-2:CENTRO(i)+2) = 2;
    tLAB(gfidx(gidx)) = 3;
    
    %plot(tLAB*100,1:500,'w')
    
    
    plot(((lidx(:,i)-1)*20)+3*61,1:500,'w')

    
    CENsig = lidx == 5;
    CENsig = 
    
    
    
    hold off

    title(num2str(probe))
  

    drawnow
    %waitforbuttonpress
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% third something level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% i think this makes the error mask for the next steps--? 
ERR = viewS - simS2;
ERRsummary = sum(sum(ERR.^2),3).^.5;
E1 = squeeze(sum(ERRsummary,1));
uE1 = mean(E1,2);
E1M = uE1 < 4;
errMask = repmat(E1M',[size(viewS,1) 1]);
errMask = repmat(errMask,[1 1 3]);
bERR = bsxfun(@times,errMask,ERR);
ERRsummary2 = squeeze(sum(sum(ERR.^2),3).^.5);
toG = [mean(ERRsummary2,1);std(ERRsummary2,1)]';
%kidx = kmeans(toG,4);
%%
close all
mag = 30;
mag = 30;
for e = 1:100:toLoad
    
    top = simC2(1,:,1,e);
    topMask = top > 0;
    
    sig1 = simC(:,:,1,e);
    sig1 = sig1 > -.83;
    
    
    sig2 = simMC(:,:,2,e);
    pl2 = sig2;
    npl2 = imfilter(pl2,fspecial('average',[300 1]),'replicate');
    npl2 = pl2 - npl2;
    npl2 = bindVec(npl2);
    
    
    
    
    
    sig3 = simC(:,:,3,e);
    pl3 = sig3;
    sig3 = sig3 < .05;
    
   
    
    sig1 = repmat(sig1,[1 size(simS,2)]);
    sig3 = repmat(sig3,[1 size(simS,2)]);
    sig1 = cat(3,sig1,sig1,sig3);
    
    
    err = rgb2gray(viewS(:,:,:,e) - simS2(:,:,:,e));
    
    img = cat(2,viewS(:,:,:,e),simS2(:,:,:,e),simMS(:,:,:,e),repmat(err,[1 1 3]));
    
    
    imshow(img,[]);
    hold on
    plot(mag*pl3,1:size(simS,1),'b');
    plot(mag*pl2+10,1:size(simS,1),'g');
    plot(mag*npl2+20,1:size(simS,1),'m');
    
    
    
    fidxc = find(imdilate(npl2,ones(101,1)) == npl2);
    mskC = npl2 > graythresh(npl2);
    kp = find(mskC(fidxc));
    fidxc = fidxc(kp);
    
    
    
    
    
    [~,fidxcn] = max(pl3);
    
    
    ZC = zeros(size(simS,1),1);
    
    
    ZCn = zeros(size(simS,1),1);
    ZC(fidxc) = 1;
    ZCn(fidxcn) = 1;
    
    plot(size(simS,2)*3*ZC,1:size(simS,1),'g');
    plot(size(simS,2)*3*ZCn,1:size(simS,1),'b');
    plot(1:numel(topMask),top*100,'c')
    drawnow
    pause(.3)
end

%% for demo
tr = 4000;
CL = {'r' 'g' 'b'};
close all
figure;
hold on
for k = 1:3
     plot(squeeze(simC(:,:,k,tr)),CL{k})
end
figure;imshow(viewS(:,:,:,tr),[]);
%% NEXT RUN - for clustering
nG = 4;
h1 = squeeze(mean(simC(1:50,:,:,:),1));
kidx = kmeans(h1',nG);
%% try different cluster method
toC = squeeze(simC);
szC = size(toC);
toC = reshape(toC,[prod(szC(1:3)) szC(4)]);
[~,CC] = PCA_FIT_FULL_T(toC,3);
gm = fitgmdist(CC',3);
kidx = gm.cluster(CC');
%%
nG = numel(unique(kidx));
for u = 1:nG
    fidx = find(kidx==u);
    fidx = fidx(randperm(numel(fidx)));
    close all
    for r = 1:10
        imshow(viewS(:,:,:,fidx(r)),[]);
        drawnow
        title([num2str(u) '--' num2str(r)]);
        waitforbuttonpress
        
    end
    
    
end
%% how many are in a group
selK = 1;
sum(kidx==selK)
subList = tmpList(kidx==selK);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first level PCA
kS = bsxfun(@times,S(:,:,:,kidx==selK),~errMask);
kS = permute(kS,[1 4 2 3]);
osz = size(kS);
kS = reshape(kS,[prod(osz(1:2)) osz(3:4)]);
%% loop over colors
clear ksimS
clear kC
clear kU
clear kE 
for e = 1:size(kS,3)
    [ksimS(:,:,e),kC(:,:,e),kU{e},kE{e}] = PCA_FIT_FULL(kS(:,:,e),1);
end
%% reshape
ksimS = reshape(ksimS,osz);
ksimS = ipermute(ksimS,[1 4 2 3]);

csz = osz;
csz(3) = 1;
ksimC = reshape(kC,csz);
ksimC = ipermute(ksimC,[1 4 2 3]);

kviewS = reshape(kS,osz);
kviewS = ipermute(kviewS,[1 4 2 3]);
kviewS = S(:,:,:,kidx==selK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first level PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
CL = {'r' 'g' 'b'};
mag = 30;
h1 = figure;
h2 = figure;
h3 = figure;
gPDF = [];
bPDF = [];
disp = false;
for e = 1:size(kviewS,4)
    img = cat(2,kviewS(:,:,:,e),ksimS(:,:,:,e));
    if disp
        figure(h1);
        imshow(img,[]);
        hold on
    end
    nsig = [];
    unsig = [];
    flippy = [1 -1 1];
    for k = 1:3
        tmpy = ksimC(:,:,k,e);
        unsig = [unsig tmpy];
        tmpy = bindVec(flippy(k)*tmpy);
        nsig = [nsig tmpy];
        tmpy = tmpy*mag;
        if disp
            plot(tmpy,1:size(ksimC,1),CL{k});
        end
    end
    
    red = squeeze(ksimC(:,:,1,e));
    red = bindVec(red);
    %{
    redMask = red > .3;
    redMask = [zeros(size(redMask)) redMask zeros(size(redMask))];
    redMask = imclearborder(redMask);
    redMask = redMask(:,2);
    redMask = [ones(size(redMask)) redMask ones(size(redMask))];
    redMask = imfill(redMask,'holes');
    redMask = redMask(:,2);
    %}
    redMask = bwlarge(red > .2);
    ridx = find(redMask);
    rsig = zeros(size(tmpy));
    rsig(ridx(1)) = 1;
    rsig(ridx(end)) = 1;
    
    %plot(rsig*size(img,2),1:size(ksimC,1),'r');
    
    
    
    
    blue = squeeze(ksimC(:,:,3,e));
    blue = blue.*redMask;
    blue = bindVec(blue);
    [NOR,bidx] = max(blue);
    bsig = zeros(size(tmpy));
    bsig(bidx) = 1;
    %plot(bsig*size(img,2),1:size(ksimC,1),'b');
    
    Rawgreen = squeeze(-ksimC(:,:,2,e));
    Rawgreen = Rawgreen.*redMask;
    green = bindVec(Rawgreen);
    gidx = find(imdilate(green,strel('disk',11)) == green);
    gsig = zeros(size(tmpy));
    gval = green(gidx);
    Rawgval = Rawgreen(gidx);
    blueNorGreen = Rawgreen * NOR^-1;
    bgval = blueNorGreen(gidx);
    %gidx = gidx(Rawgval > .4);
    %gidx = gidx(gval > .5 & Rawgval > .4);
    %gidx = gidx(gval > .5 & Rawgval > .85);
    %gidx = gidx(gval > .5 & Rawgval > .85);
    gidx = gidx(gval > .5 & bgval >2* 1/6);
    
    
    
    gsig(gidx) = 1;
    
    if bidx > 250
        gsig = flipdim(gsig,1);
        bsig = flipdim(bsig,1);
    end
    
    gPDF = [gPDF gsig];
    bPDF = [bPDF bsig];
    
    %plot(gsig*size(img,2),1:size(ksimC,1),'g');
    
    %{
    figure(h2);
    plot(nsig);
    
    
    figure(h3);
    plot(unsig);
    
    figure(h1);
    hold off
    %}
    if disp
        drawnow
        waitforbuttonpress
    end
    e
    N(e) = sum(gsig);
end


%%
close all
plot(mean(gPDF,2))
%%
n = sum(gPDF,1);
close all
totC = 4;
UF = [];
UFB = [];
for e = 1:totC
    
    
    fidx = find(n == e);
    tmp = gPDF(:,fidx);
    UF(:,e) = mean(tmp,2);
    UF(:,e) = imfilter(UF(:,e),fspecial('average',[31 1]));
    UF(:,e) = UF(:,e) * sum(UF(:,e))^-1;
    
    tmp = bPDF(:,fidx);
    UFB(:,e) = mean(bPDF,2);
    UFB(:,e) = imfilter(UFB(:,e),fspecial('average',[31 1]));
    UFB(:,e) = UFB(:,e) * sum(UFB(:,e))^-1;
    
    
end
for e = 1:totC
    figure;
    plot(UF(:,e))
    hold on
    plot(UFB(:,e))
end
%%
outTable = table;
outTable.name = subList';
for e = 1:size(gPDF,1)
    outTable.(['P' num2str(e)]) = gPDF(e,:)';
    e
end
%%
writetable(outTable,'/mnt/tetra/nate/try2.csv')
%%
for e = 1:numel(H)
    rm(e) = contains(lower(H{e}),'bounding');
end
H(find(rm)) = [];
D(:,find(rm)) = [];
%%

background = hmm_node('background');
default = hmm_node('default');
centromere = hmm_node('centromere');
foci = hmm_node('foci');



hmm = my_hmm();
hmm.addNode(background);
hmm.addNode(default);
hmm.addNode(centromere);
hmm.addNode(foci);
% FIX
hmm.dn = ones(numel(wC),1);



backgroundCov = cov(wC{m}(:,preIdx)');
backgroundU{1,m} = mean(wC{m}(:,preIdx)');
    
    
    tranCov{2,m} = cov(wC{m}(:,tranIdx)');
    tranU{2,m} = mean(wC{m}(:,tranIdx)');
    
    postCov{3,m} = cov(wC{m}(:,postIdx)');
    postU{3,m} = mean(wC{m}(:,postIdx)');




















