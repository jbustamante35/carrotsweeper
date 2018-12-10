inFilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scan for new images
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(inFilePath,FileList,FileExt,verbose);
%% run
basePath = '//mnt/snapper/nate/myDev/maizeWhole/';
mkdir(basePath);
for e = 1:numel(SET)
    % generate save name
    [pth nm ext] = fileparts(SET{e}{1});
    pth = strrep(pth,filesep,'-');
    saveName = [pth '.mat'];
    
    tic
    cropMaizeImage(SET{e},saveName,basePath);   
    toc
end
%% iterative load the whole image
PER = 1;
for e = 1:10
    
inFilePath = '/mnt/spaldingdata/nate/mirror_images/beta_miniMaize/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = gdig(inFilePath,FileList,FileExt,verbose);



Itmp = imread(SET{1});
Itmp = imresize(Itmp,PER);
Itmp = Itmp(round(end-end*.75):round(end-end*.25),:);
sz = size(Itmp);
M = zeros(numel(SET),numel(Itmp(:)));
rm = [];



for i = 1:numel(SET)
    try
        I = imread(SET{i});
        I = imresize(I,PER);
        I = I(round(end-end*.75):round(end-end*.25),:);
        M(i,:) = I(:);
        i
        numel(SET)
    catch
        rm = [rm i];
    end
end


M(rm,:) = [];
M = double(M);
[S C U E L ERR LAM] = PCA_FIT_FULL(M,200);

close all
for i = 1:size(E,2)
    
    imshow(reshape(E(:,i),sz),[]);
    drawnow
    pause(1)
end

for i = 1:size(S,1)    
    imshow(reshape(S(i,:),sz),[]);
    drawnow
    pause(1)
end


pause(1000)
end
%% load the second half via sliding windows
PER = 1;
for e = 1:10

    inFilePath = '/mnt/spaldingdata/nate/mirror_images/beta_miniMaize/';
    FileList = {};
    FileExt = {'tiff','TIF','tif'};
    verbose = 1;
    SET = gdig(inFilePath,FileList,FileExt,verbose);


   
    Itmp = imread(SET{1});
    Itmp = imresize(Itmp,PER);
    Itmp = Itmp(round(end-end*.75):round(end-end*.25),round(end/2):end);
    Itmp = im2col(Itmp,[51 51],'sliding');
    sz = size(Itmp);
    M = zeros(sz(1),numel(SET)*sz(2));
    rm = [];



for i = 1:numel(SET)
    try
        I = imread(SET{i});
        I = imresize(I,PER);
        I = I(round(end-end*.75):round(end-end*.25),round(end/2):end);
        I = im2col(I,[51 51],'sliding');
        M(,:) = I(:);
        i
        numel(SET)
    catch
        rm = [rm i];
    end
end


M(rm,:) = [];
M = double(M);
[S C U E L ERR LAM] = PCA_FIT_FULL(M,200);

close all
for i = 1:size(E,2)
    
    imshow(reshape(E(:,i),sz),[]);
    drawnow
    pause(1)
end

for i = 1:size(S,1)    
    imshow(reshape(S(i,:),sz),[]);
    drawnow
    pause(1)
end


pause(1000)
end
%%  load the whole along the third
   
inFilePath = '/mnt/spaldingdata/nate/mirror_images/beta_miniMaize/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = gdig(inFilePath,FileList,FileExt,verbose);


PER = 1;
Itmp = imread(SET{1});
Itmp = imresize(Itmp,PER);
Itmp = Itmp(round(end-end*.75):round(end-end*.25),:);
sz = size(Itmp);
M = zeros([sz numel(SET)]);
rm = [];


rm = [];
for i = 1:numel(SET)
    try
        I = imread(SET{i});
        I = imresize(I,PER);
        I = I(round(end-end*.75):round(end-end*.25),:);
        M(:,:,i) = I;
        i
        numel(SET)
    catch
        rm = [rm i];
    end
end
M(:,:,rm) = [];
%%  dig for mat data   
inFilePath = '//mnt/snapper/nate/myDev/maizeFirstFrame/';
FileList = {};
FileExt = {'mat'};
verbose = 1;
SET = gdig(inFilePath,FileList,FileExt,verbose);
%% load better
MIP = [];
NM = [];
CB = [];
SEG = [];
EE = [];
DV = [];
for e = 1:500%numel(SET)
    data = load(SET{e});
    %{
    for tm = 1:size(data.tmpCurve.S,3)
        imshow(data.tmpCurve.S(:,:,tm),[]);
        drawnow
    end
    %}
    
    MIP = cat(3,MIP,data.tmpCurve.S);
    EE = cat(3,EE,data.tmpCurve.E);
    NM = [NM;e*ones(size(data.tmpCurve.S,3),1)];
    CB = [CB;data.tmpCurve.data'];
    DV = [DV;bsxfun(@minus,data.tmpCurve.data,mean(data.tmpCurve.data,2))'];
    SEG = [SEG;reshape(data.tmpCurve.segs,[size(data.tmpCurve.segs,1)*size(data.tmpCurve.segs,2) size(data.tmpCurve.segs,3)])'];
    %SEG = cat(3,SEG,data.tmpCurve.segs);
    e
end
MIP = reshape(MIP,[size(MIP,1)*size(MIP,2) size(MIP,3)]);
%% PCA on data
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(MIP',10);
%% PCA on data
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(SEG,10);
%% display tip
close all
UQ = unique(NM);
h1 = figure;
h2 = figure;
h3 = figure;
for u = 1:numel(UQ)
    figure(h1);    
    fidx = find(NM==UQ(u));
    Psig = pC(fidx,1:3);    
    plot(Psig)
    title('Patch PCA');
    drawnow
    
    figure(h3);    
    fidx = find(NM==UQ(u));
    Ssig = sC(fidx,1:3);
    plot(Ssig)
    title('Segment PCA');
    drawnow
    
    figure(h2);    
    tmpC = CB(fidx,:);
    plot(tmpC(:,1),tmpC(:,2));
    
    E = EE(:,:,fidx);
    
    hold on
    [~,sidx1] = max(Psig(:,2));
    [~,sidx2] = min(Ssig(:,1));
    plot(tmpC(sidx1,1),tmpC(sidx1,2),'ro');
    plot(tmpC(sidx2,1),tmpC(sidx2,2),'g*');
    quiver(tmpC(:,1),tmpC(:,2),squeeze(E(1,2,:)),squeeze(E(2,2,:)));
    axis equal
    hold off
    waitforbuttonpress;
end
%% segment on patches via gmm
k = 4;
ndims = 2;
opt = statset('display','iter','MaxIter',500);
obj = gmdistribution.fit(pC(:,1:ndims),k,'Options',opt);
%% segment via kmeans
k = 4;
kidx = kmeans(pC(:,1:4),k);
%% segment on segments
k = 3;
opt = statset('display','iter');
obj = gmdistribution.fit(sC(:,1:4),k,'Options',opt);
%% display segmentation from gmm
close all
UQ = unique(NM);
h1 = figure;
CL = {'r' 'g' 'b' 'k' 'y' 'm'};
for u = 1:numel(UQ)
    figure(h1);    
    fidx = find(NM==UQ(u));
    Psig = pC(fidx,1:ndims);    
    
    
    y = [];
    for g = 1:size(obj.mu,1)
        y(:,g) = mvnpdf(Psig,obj.mu(g,:),obj.Sigma(:,:,g));
    end
    [~, gidx] = max(y,[],2);
    
    
    
    
    figure(h2);    
    tmpC = CB(fidx,:);
    plot(tmpC(:,1),tmpC(:,2));
    % get tangent and normal space
    E = EE(:,:,fidx);    
    hold on
    
    %[~,sidx1] = max(Psig(:,2));    
    %plot(tmpC(sidx1,1),tmpC(sidx1,2),'b*');    
    quiver(tmpC(:,1),tmpC(:,2),squeeze(E(1,2,:)),squeeze(E(2,2,:)));
    for g = 1:size(obj.mu,1)
        fidx = find(gidx==g);
        plot(tmpC(fidx,1),tmpC(fidx,2),[CL{g} 'o']);
    end
    axis equal
    hold off
    waitforbuttonpress;
end
%% segment from gmm and ...
y = [];
for g = 1:size(obj.mu,1)
    y(:,g) = mvnpdf(pC(:,1:ndims),obj.mu(g,:),obj.Sigma(:,:,g));
end
[~, gidx_master_Y] = max(y,[],2);
%% get parameters for distributions
tipGrp = 4;
rootGrps = 1;
transGrps = 3;
kernelGrp = 2;
ndims = 2;
tip_mean_patch = mean(pC(gidx_master_Y==tipGrp,1:ndims));
tip_cov_patch = cov(pC(gidx_master_Y==tipGrp,1:ndims));
root_mean_patch = mean(pC(gidx_master_Y==rootGrps,1:ndims));
root_cov_patch = cov(pC(gidx_master_Y==rootGrps,1:ndims));
trans_mean_patch = mean(pC(gidx_master_Y==transGrps,1:ndims));
trans_cov_patch = cov(pC(gidx_master_Y==transGrps,1:ndims));
kernel_mean_patch = mean(pC(gidx_master_Y==kernelGrp,1:ndims));
kernel_cov_patch = cov(pC(gidx_master_Y==kernelGrp,1:ndims));

tip_mean_cord = mean(DV(gidx_master_Y==tipGrp,:));
tip_cov_cord = cov(DV(gidx_master_Y==tipGrp,:));
root_mean_cord = mean(DV(gidx_master_Y==rootGrps,:));
root_cov_cord = cov(DV(gidx_master_Y==rootGrps,:));
trans_mean_cord = mean(DV(gidx_master_Y==transGrps,:));
trans_cov_cord = cov(DV(gidx_master_Y==transGrps,:));
kernel_mean_cord = mean(DV(gidx_master_Y==kernelGrp,:));
kernel_cov_cord = cov(DV(gidx_master_Y==kernelGrp,:));



%% generate HMM
tip_node = hmm_node('tip');
root_node_upper = hmm_node('root_upper');
t1_node = hmm_node('ut');
kernel_node = hmm_node('kernel_node');
t2_node = hmm_node('lt');
root_node_lower = hmm_node('root_lower');
kernel_node_end = hmm_node('e_kerel_node');



tf1 = heavisideTransitionFunction(20,@(x,y)lt(x,y));
tf2 = heavisideTransitionFunction(20,@(x,y)ge(x,y));
kernel_to_transZone = constantTransitionFunction(.1);
kernel_to_kernel = constantTransitionFunction(.9);
kernelEnd_to_kernelEnd = constantTransitionFunction(1);
root_to_root = constantTransitionFunction(.75);
root_to_outofToot = constantTransitionFunction(.25);
tip_to_tip = constantTransitionFunction(.75);
tip_to_lower_root = constantTransitionFunction(.25);

kernel_node.attachNode(kernel_node,kernel_to_kernel);
kernel_node.attachNode(t1_node,kernel_to_transZone);

t1_node.attachNode(t1_node,tf1);
t1_node.attachNode(root_node_upper,tf2);

root_node_upper.attachNode(root_node_upper,root_to_root);
root_node_upper.attachNode(tip_node,root_to_outofToot);

%tip_node.attachNode(tip_node,tf1);
%tip_node.attachNode(root_node_lower,tf2);

tip_node.attachNode(tip_node,tip_to_tip);
tip_node.attachNode(root_node_lower,tip_to_lower_root);


root_node_lower.attachNode(root_node_lower,root_to_root);
root_node_lower.attachNode(t2_node,root_to_outofToot);

t2_node.attachNode(t2_node,tf1);
t2_node.attachNode(kernel_node_end,tf2);

kernel_node_end.attachNode(kernel_node_end,kernelEnd_to_kernelEnd);






% attach distributions for patches
tipD = myProb(tip_mean_patch,tip_cov_patch);
tip_node.attachDistribution(tipD,1);
rootD = myProb(root_mean_patch,root_cov_patch);
root_node_lower.attachDistribution(rootD,1);
root_node_upper.attachDistribution(rootD,1);
transD = myProb(trans_mean_patch,trans_cov_patch);
t1_node.attachDistribution(transD,1);
t2_node.attachDistribution(transD,1);
kernelD = myProb(kernel_mean_patch,kernel_cov_patch);
kernel_node.attachDistribution(kernelD,1);
kernel_node_end.attachDistribution(kernelD,1);
% attach distributions for cordinates
tipD_cord = myProb(tip_mean_cord,tip_cov_cord);
tip_node.attachDistribution(tipD_cord,2);
rootD_cord = myProb(root_mean_cord,root_cov_cord);
root_node_lower.attachDistribution(rootD_cord,2);
root_node_upper.attachDistribution(rootD_cord,2);
transD_cord = myProb(trans_mean_cord,trans_cov_cord);
t1_node.attachDistribution(transD_cord,2);
t2_node.attachDistribution(transD_cord,2);
kernelD_cord = myProb(kernel_mean_cord,kernel_cov_cord);
kernel_node.attachDistribution(kernelD_cord,2);
kernel_node_end.attachDistribution(kernelD_cord,2);
% create hmm
hmm = my_hmm();
hmm.addNode(tip_node);
hmm.addNode(root_node_upper);
hmm.addNode(t1_node);
hmm.addNode(kernel_node);
hmm.addNode(t2_node);
hmm.addNode(root_node_lower);
hmm.addNode(kernel_node_end);
hmm.dn = [1 1];
tm = hmm.buildTransitionMatrix([]);
%% check space distribution map
[n1 n2] = ndgrid(-300:300,-300:300);
N = [n1(:) n2(:)];
PDF_root_cord = rootD_cord.getProb(N',1);
PDF_root_cord = reshape(PDF_root_cord,size(n1));
PDF_trans_cord = transD_cord.getProb(N',1);
PDF_trans_cord = reshape(PDF_trans_cord,size(n1));
PDF_kernel_cord = kernelD_cord.getProb(N',1);
PDF_kernel_cord = reshape(PDF_kernel_cord,size(n1));
PDF_tip_cord = tipD_cord.getProb(N',1);
PDF_tip_cord = reshape(PDF_tip_cord,size(n1));
close all
imshow(cat(3,bindVec(PDF_root_cord)+bindVec(PDF_tip_cord),bindVec(PDF_trans_cord)+bindVec(PDF_tip_cord),bindVec(PDF_kernel_cord)+bindVec(PDF_tip_cord)),[])
figure;
hold on
plot(tip_mean_cord(1),tip_mean_cord(2),'b*')
plot(root_mean_cord(1),root_mean_cord(2),'g*')
plot(trans_mean_cord(1),trans_mean_cord(2),'k*')
plot(kernel_mean_cord(1),kernel_mean_cord(2),'r*')

%% clear class hmm and node_hnm
clear tip_node root_node_upper t1_node kernel_node t2_node root_node_lower kernel_node_end hmm kernel_to_transZone kernel_to_kernel kernelEnd_to_kernelEnd root_to_root root_to_outofToot
clear tf1 tf2
clear class my_hmm 
clear class hmm_node
clear class myProb
clear class myTransitionFunction
clear class constantTransitionFunction
clear class heavisideTransitionFunction
%% display segmentation from hmm
close all
UQ = unique(NM);
h1 = figure;
CL = {'r' 'b' 'g' 'k' 'y' 'c' 'm' 'c'};
recalc = 0;
for u = 1:numel(UQ)
    figure(h1);    
    fidx = find(NM==UQ(u));
    Psig_patch = pC(fidx,1:ndims);    
    Psig_cord = DV(fidx,1:ndims);
    Psig_tot = [Psig_patch';Psig_cord'];
    observation_labels = [ones(size(Psig_patch,2),1);2*ones(size(Psig_cord,2),1)];
    
    
    if recalc
        gidx{u} = hmm.Viterbi(Psig_tot,observation_labels,4);
    end
       
    tmpC = CB(fidx,:);
    plot(tmpC(:,1),tmpC(:,2));    
    % get tangent and normal space
    E = EE(:,:,fidx);    
    hold on
    
    
    quiver(tmpC(:,1),tmpC(:,2),squeeze(E(1,2,:)),squeeze(E(2,2,:)),'b');
    UQg = unique(gidx{u});
    for g = 1:numel(UQg)
        fidx = find(gidx{u}==UQg(g));
        plot(tmpC(fidx,1),tmpC(fidx,2),[CL{UQg(g)} '*']);
    end
    %plot(tmpC(1,1),tmpC(1,2),'ro');
    %plot(tmpC(2,1),tmpC(2,2),'go');
    axis equal
    hold off
    drawnow
    waitforbuttonpress;
end
%% generate update information
close all
UQ = unique(NM);
UD = {};
for u = 1:numel(UQ)    
    fidx = find(NM==UQ(u));
    Psig_patch = pC(fidx,1:ndims);    
    Psig_cord = DV(fidx,1:ndims);
    Psig_tot = [Psig_patch';Psig_cord'];
    observation_labels = [ones(size(Psig_patch,2),1);2*ones(size(Psig_cord,2),1)];
    gidx = hmm.Viterbi(Psig_tot,observation_labels,4);
    UD{u} = [gidx;Psig_tot];
    u
end
%% remove offending data
[idx] = hmm.detectIllegalTransitions(UD);
%% update
observation_labels = [ones(ndims,1);2*ones(size(DV,2),1)];
hmm.update(UD,observation_labels);
%% simulate
[n1 n2] = ndgrid(-300:300,-300:300);
N = [n1(:) n2(:)];
mPDF = zeros(size(n1));
for s = 1:numel(hmm.NodeList)
    PDF = hmm.NodeList{s}.Distribution{2}.getProb(N',1);
    PDF = reshape(PDF,size(n1));
    PDF = bindVec(PDF);
    imshow(PDF,[]);
    title(hmm.NodeList{s}.StateName);
    waitforbuttonpress
    mPDF = mPDF + PDF;
end
%%
imshow(mPDF,[]);
%% loop
for e = 1:100

    e
    
    
    close all
    UQ = unique(NM);
    UD = {};
    for u = 1:numel(UQ)    
        fidx = find(NM==UQ(u));
        Psig_patch = pC(fidx,1:ndims);    
        Psig_cord = DV(fidx,1:ndims);
        Psig_tot = [Psig_patch';Psig_cord'];
        observation_labels = [ones(size(Psig_patch,2),1);2*ones(size(Psig_cord,2),1)];
        gidx = hmm.Viterbi(Psig_tot,observation_labels,4);
        UD{u} = [gidx;Psig_tot];
        u
    end
    
    
    
    observation_labels = [ones(ndims,1);2*ones(size(DV,2),1)];
    hmm.update(UD,observation_labels);
    
end