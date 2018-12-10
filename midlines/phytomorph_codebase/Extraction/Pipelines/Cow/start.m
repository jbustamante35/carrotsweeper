FilePath = '/mnt/tetra/nate/cow/tif2/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% scan view
close all
for e = 1:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow
end
%%
close all
cnt = 1;
SKIP = 30;
loadVec = 1:SKIP:numel(FileList);
loadVec = 1:900;
S = imread(FileList{e});
SNIP = 30;
S((end-SNIP):end,:,:) = [];
per = .45;
M = imresize(S,per);
g = rgb2gray(M);
S = zeros([size(S) numel(loadVec)]);
G = zeros([size(S,1) size(S,2) numel(loadVec)]);
g = zeros([size(g,1),size(g,2),numel(loadVec)]);
M = zeros([size(M) numel(loadVec)]);
lList = FileList(loadVec);



parfor e = 1:numel(lList)
    tmp = double(imread(lList{e}))/255;
    tmp((end-SNIP):end,:,:) = [];
    
    tmp = imhistmatch(tmp,HM);
    %tmp = imfilter(tmp,fspecial('gaussian',[21 21],4),'replicate');
    S(:,:,:,e) = tmp;
    G(:,:,e) = rgb2gray(S(:,:,:,e));
    M(:,:,:,e) = imresize(S(:,:,:,e),per);
    g(:,:,e) = rgb2gray(M(:,:,:,e));
    %{
    imshow(S(:,:,:,e),[]);
    title(num2str(e))
    drawnow
    %}
    e
end
%%
HM = mean(S,4);
%% decompose gray
mSZ = size(g);
vec = reshape(g,[prod(mSZ(1:2)) mSZ(3)]);
[vec,mu,std] = zscore(vec,1,2);
[U,E,L] = PCA_FIT_FULL_Tws(vec,5);
%% decompose color
mSZ = size(M);
vec = reshape(M,[prod(mSZ(1:3)) mSZ(4)]);
[U,E,L] = PCA_FIT_FULL_Tws(vec,2);
%% decompose color channels
clear U E L
for k = 1:3
    vec = squeeze(M(:,:,k,:));
    mSZ = size(vec);
    vec = reshape(vec,[prod(mSZ(1:2)) mSZ(3)]);
    [U{k},E{k},L{k}] = PCA_FIT_FULL_Tws(vec,5);
end
%% upscale color
vU = reshape(U,mSZ(1:3));
imshow(vU,[])
vE = reshape(E(:,1),mSZ(1:3));
imshow(bindVec(vE),[])
%% upscale gray
close all
sSZ = size(G);
vU = reshape(U,mSZ(1:2));
vU = imresize(vU,sSZ(1:2));
imshow(vU,[]);
waitforbuttonpress
nV = vU(:);
nE = [];
for e = 1:5
    vE = reshape(E(:,e),mSZ(1:2));
    imshow(vE,[]);
    drawnow
    waitforbuttonpress
    vE = imresize(vE,sSZ(1:2));
    imshow(vE,[]);
    drawnow
    waitforbuttonpress
    nE(:,e) = vE(:);
end
%% upscale color 2
clear nV nE
for k = 1:3
    close all
    mSZ = size(g);
    sSZ = size(G);
    vU = reshape(U{k},mSZ(1:2));
    vU = imresize(vU,sSZ(1:2));
    imshow(vU,[]);
    %waitforbuttonpress
    nV{k} = vU(:);
    for e = 1:2
        vE = reshape(E{k}(:,e),mSZ(1:2));
        imshow(vE,[]);
        vE = imresize(vE,sSZ(1:2));
        imshow(vE,[]);
        drawnow
        %waitforbuttonpress
        nE{k}(:,e) = vE(:);
    end
end
close all
%% mu std
mSZ = size(g);
vec = reshape(g,[prod(mSZ(1:2)) mSZ(3)]);
[vec,mu,std] = zscore(vec,1,2);
mu = reshape(mu,mSZ(1:2));
mu = imresize(mu,sSZ(1:2));
imshow(mu,[]);
waitforbuttonpress
mu = mu(:);
std = reshape(std,mSZ(1:2));
std = imresize(std,sSZ(1:2));
imshow(std,[])
waitforbuttonpress
std = std(:);
%% upscale color
sSZ = size(S);
vU = reshape(U,mSZ(1:3));
vU = imresize(vU,sSZ(1:2));
nV = vU(:);
for e = 1
    vE = reshape(E(:,e),mSZ(1:3));
    vE = imresize(vE,sSZ(1:2));
    nE(:,e) = vE(:);
end
%%
close all
C = [];

close all
SKIP = 50;
loadVec = 1:SKIP:numel(FileList);
lList = FileList(loadVec);

for e = 1:size(S,4)
    
    tmp = double(imread(lList{e}))/255;
    tmp((end-SNIP):end,:,:) = [];
    tmp = imhistmatch(tmp,HM);
    %{
    tmp = G(:,:,e);
    C(e,:) = PCA_REPROJ_T(tmp(:),nV,nE);
    SIM = PCA_BKPROJ_T(C(e,:)',nE,nV);
    SIM = reshape(SIM,size(tmp));
    %}
    %tmp = S(:,:,:,e);
    tSIM = [];
    for k = 1:3
        vec = squeeze(tmp(:,:,k));
        mSZ = size(vec);
        C = PCA_REPROJ_T(vec(:),nE{k},nV{k});
        SIM = PCA_BKPROJ_T(C,nE{k},nV{k});
        tSIM = [tSIM SIM];
    end
    tSIM = reshape(tSIM,size(tmp));
    delta = tmp - tSIM;
    delta = sum(delta.*delta,3).^.5;
    mask = abs(delta) > .15;
    mask = bwareaopen(mask,20);
    v = bsxfun(@times,tmp,mask);
    %{
    imshow(v,[]);
    %}
    %{
    C = PCA_REPROJ_T(tmp(:),nV,nE);
    C = PCA_BKPROJ_T(C,nV,nE);
    C = reshape(C,size(tmp));
    %}
    %{
    for k = 1:size(C,3)
        C(:,:,k) = bindVec(C(:,:,k));
    end
    %}
    imshow(delta,[]);
    %imshow([tmp tSIM],[]);
    drawnow
end
%% build up delta stats
close all
C = [];

close all
SKIP = 50;
loadVec = 1:SKIP:numel(FileList);
lList = FileList(loadVec);

toG = 900;
gu = [];
su = [];
%disD = [];
clear std;
for e = 1:size(S,4)
    
    tmp = double(imread(lList{e}))/255;
    tmp((end-SNIP):end,:,:) = [];
    tmp = imhistmatch(tmp,HM);
    %{
    tmp = G(:,:,e);
    C(e,:) = PCA_REPROJ_T(tmp(:),nV,nE);
    SIM = PCA_BKPROJ_T(C(e,:)',nE,nV);
    SIM = reshape(SIM,size(tmp));
    %}
    %tmp = S(:,:,:,e);
    tSIM = [];
    for k = 1:3
        vec = squeeze(tmp(:,:,k));
        mSZ = size(vec);
        C = PCA_REPROJ_T(vec(:),nE{k},nV{k});
        SIM = PCA_BKPROJ_T(C,nE{k},nV{k});
        tSIM = [tSIM SIM];
    end
    
    
    
    tSIM = reshape(tSIM,size(tmp));
    deltaC = tmp - tSIM;
    delta = sum(deltaC.*deltaC,3).^.5;
    
    if loadVec(e) < toG
        disD = cat(3,disD,delta);
        gu = mean(disD,3);
        su = std(disD,1,3);
        %{
        gu = mean(cat(3,gu,delta));
        gu = mean([gu mean(delta(:))]);
        su = mean([su std(delta(:),1,1)]);
        %}
    end
    
    z = (delta - gu).*su.^-1;
    
    mask = abs(delta) > .15;
    mask = abs(z) > 2.1;
    mask = bwareaopen(mask,50);
    v = bsxfun(@times,tmp,mask);
    %{
    imshow(v,[]);
    %}
    %{
    C = PCA_REPROJ_T(tmp(:),nV,nE);
    C = PCA_BKPROJ_T(C,nV,nE);
    C = reshape(C,size(tmp));
    %}
    %{
    for k = 1:size(C,3)
        C(:,:,k) = bindVec(C(:,:,k));
    end
    %}
    imshow(v,[]);
    title(num2str(loadVec(e)))
    %imshow([tmp tSIM],[]);
    drawnow
end
%% page load and apply
close all
SKIP = 50;
loadVec = 1:SKIP:numel(FileList);
lList = FileList(loadVec);
h1 = figure;
h2 = figure;
for e = 1:numel(lList)
    tmp = double(imread(lList{e}))/255;
    tmp((end-SNIP):end,:,:) = [];
    %ftmp = imfilter(tmp,fspecial('gaussian',[21 21],3),'replicate');
   
    gtmp = rgb2gray(ftmp);
    szGG = size(gtmp);
    gtmp = gtmp(:);
    gtmp = (gtmp - mu).*std.^-1;
    C = PCA_REPROJ_T(gtmp(:),nE,nV);
    SIM = PCA_BKPROJ_T(C,nE,nV);
    SIM = reshape(SIM,szGG);
    
    mask = abs([reshape(gtmp,szGG)-SIM]) > 15;
    ctmp = bsxfun(@times,mask,tmp);
    imshow(ctmp,[]);
    drawnow
    
    %{
    figure(h1);
    imshow([SIM tmp tmp-SIM],[]);
    title(num2str(lList{e}))
    drawnow
    %}
    %{
    figure(h2)
    imshow(reshape(gtmp,szGG)-SIM,[-20 20]);
    waitforbuttonpress
    %}
    %{
    
    imshow(S(:,:,:,e),[]);
    title(num2str(e))
    drawnow
    %}
    e
end