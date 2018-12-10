FilePath = '/mnt/snapper/nate/forRichard/phytoMorph_Image_Phenomics_Tool_Kit_analysis1-2017-07-10-21-23-22.7/';
FilePath = '/mnt/snapper/nate/forRichard/phytoMorph_Image_Phenomics_Tool_Kit_analysis1-2017-07-13-17-57-01.6/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
COT = readtext('/home/nate/Downloads/cots.csv');
cD = cell2mat(COT(2:321,:))';
H = COT(1,:);



ncD = bsxfun(@minus,cD,mean(cD(:,1:30),2));


%% 1: load data
cnt = 1;
data = [];
N = 320;
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'rootTipSignal'))
        tmp = csvread(FileList{e});
        data(cnt).d = tmp(1:N)';
        [pth,nm,ext] = fileparts(FileList{e});
        data(cnt).n = nm;
        tmp = csvread(strrep(FileList{e},'rootTipSignal','kymogrpah'));
        data(cnt).k = tmp(:,1:N);
        tmp = csvread(strrep(FileList{e},'rootTipSignal','growthRate'));
        data(cnt).g = tmp(1:N)';
        cnt = cnt + 1;
    end
end
% 2: normalize data
close all
D = [data.d]';
G = [data.g]';
K = [];
for e = 1:numel(data)
    K = cat(3,K,data(e).k);
end
init = mean(D(:,1:30),2);
nD = bsxfun(@minus,D,init);
plot(nD');
rmidx = find(abs(nD(:,1)) > .1);



data(rmidx) = [];
nD(rmidx,:) = [];
G(rmidx,:) = [];
D(rmidx,:) = [];
init = mean(G(:,1:3),2);
nG = bsxfun(@minus,G,init);
%nG = nG(:,50:end);
nG(:,1:50) = 0;
nG = bsxfun(@minus,nG,nG(:,1));
nG = imfilter(nG,fspecial('average',[1 20]),'replicate');

plot(nD');
figure;plot(nG');
%{
close all

for e = 1:size(nG,1)
    plot(nG(e,:))
    drawnow
    pause(.2)
end

sG = nG(:,100:200);

for e = 1:size(sG,1)
    plot(sG(e,:))
    drawnow
    pause(.2)
end
%}
%% 3: fitdata
[S C U E L ERR LAM] = PCA_FIT_FULL(nD,5);
uuD = [ncD ;nD];
j1 = mean(ncD,1);
j2 = mean(nD,1);
uuD = [bsxfun(@minus,ncD,j1);bsxfun(@minus,nD,j2)];

[uS uC uU uE uL uERR uLAM] = PCA_FIT_FULL(uuD,4);

[uS1 uC1 uU1 uE1 uL1 uERR1 uLAM1] = PCA_FIT_FULL(bsxfun(@minus,ncD,j1),5);
[uS2 uC2 uU2 uE2 uL2 uERR2 uLAM2] = PCA_FIT_FULL(bsxfun(@minus,nD,j2),5);
su = sign(diag(uE1'*uE2));
for c = 1:size(uC2,2)
    uC2(:,c) = su(c)*uC2(:,c);
    uE2(:,c) = su(c)*uE2(:,c);
end
uC = [uC1;uC2];
[nnS nnC nnU nnE nnL nnERR nnLAM] = PCA_FIT_FULL(D,5);
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL(nG,2);
sz = size(K);
rK = reshape(K,[prod(sz(1:2)) sz(3)]);
[kS kC kU kE kL kERR kLAM] = PCA_FIT_FULL(rK,5);
kS = reshape(kS,sz);

%% direct compare
close all
midx = [];
widx = [];
for e = 1:numel(data)
    tmpName = lower(data(e).n);
    if  ~isempty(strfind(tmpName,'gc')) & isempty(strfind(tmpName,'etio')) & isempty(strfind(tmpName,'chitin')) & (~isempty(strfind(tmpName,'4.11')) | ~isempty(strfind(tmpName,'230')) | ~isempty(strfind(tmpName,'231')))
        midx = [midx e];
        tmpName
    end
    if  ~isempty(strfind(tmpName,'gc')) & isempty(strfind(tmpName,'etio')) & isempty(strfind(tmpName,'chitin')) & ~(~isempty(strfind(tmpName,'4.11')) | ~isempty(strfind(tmpName,'230')) | ~isempty(strfind(tmpName,'231')))
        widx = [widx e];
        tmpName
    end
    %tmpName
end
SIG = S;
%SIG = gS;
%SIG = nnS;
SIG = nD;

mU = mean(SIG(midx,:));
mS = std(SIG(midx,:),1,1)*numel(midx)^-.5;
wU = mean(SIG(widx,:));
wS = std(SIG(widx,:),1,1)*numel(widx)^-.5;
close all
errorbar(wU,wS,'k');
hold on
errorbar(mU,mS,'r');
LEG{1} = 'wildType';
LEG{2} = 'mut';
legend(LEG);

plot(SIG(widx,:)','k--')

plot(SIG(midx,:)','r--')
figure;
plot(SIG(widx,:)','k--')
hold on
plot(SIG(midx,:)','r--')
%% root only

G = [zeros(numel(widx),1);ones(numel(midx),1)];
DD = [uC(widx+18,:);uC(midx+18,:)];

%{
[lambda] = myLDA(DD,G);
[lambda] = mynLDA(DD,G,1,2);
%}
h = [];
p = [];
for l = 1:size(lambda,2)
    wv = uC(widx+18,:)*lambda(:,l);
    mv = uC(midx+18,:)*lambda(:,l);
    [h(l) p(l)] = ttest2(wv,mv);
    
    figure;
    ksdensity(wv);
    hold all
    ksdensity(mv);
end

splitVec = PCA_BKPROJ(lambda',uE,uU);


SIG = uS;
%SIG = gS;
%SIG = nnS;
%SIG = D;

mU = mean(SIG(midx+18,:));
mS = std(SIG(midx+18,:),1,1)*numel(midx)^-.5;
wU = mean(SIG(widx+18,:));
wS = std(SIG(widx+18,:),1,1)*numel(widx)^-.5;
figure;
errorbar(wU,wS,'k');
hold on
errorbar(mU,mS,'r');
for l = 1:size(splitVec,1)
    plot(2*(splitVec(l,:)-U),'b');
    
end
LEG{1} = 'wildType';
LEG{2} = 'mut';
legend(LEG);

plot(SIG(widx,:)','k--')

plot(SIG(midx,:)','r--')
figure;
plot(SIG(widx,:)','k--')
hold on
plot(SIG(midx,:)','r--')
%% cots only

G = [zeros(9,1);ones(9,1)];
DD = [uC(1:9,:);uC(10:18,:)];

%{
[lambda] = myLDA(DD,G);
[lambda] = mynLDA(DD,G,1,5);
%}
h = [];
p = [];
for l = 1:size(lambda,2)
    wv = uC(1:9,:)*lambda(:,l);
    mv = uC(10:18,:)*lambda(:,l);
    [h(l) p(l)] = ttest2(wv,mv);
end

close all
ksdensity(wv);
hold all
ksdensity(mv);
splitVec = PCA_BKPROJ(lambda',uE,uU);


SIG = uS;
%SIG = gS;
%SIG = nnS;
%SIG = D;

mU = mean(SIG(10:18,:));
mS = std(SIG(10:18,:),1,1)*numel(10:18)^-.5;
wU = mean(SIG(1:9,:));
wS = std(SIG(1:9,:),1,1)*numel(1:9)^-.5;
close all
errorbar(wU,wS,'k');
hold on
errorbar(mU,mS,'r');
for l = 1:size(splitVec,1)
    plot(2*(splitVec(l,:)-U),'b');
    
end
LEG{1} = 'wildType';
LEG{2} = 'mut';
legend(LEG);


%% ind tissue


G = [zeros(numel(widx)+9,1);ones(numel(midx)+9,1)];
DD = [uC(1:9,:);uC(widx+18,:);uC(10:18,:);uC(midx+18,:)];


[lambda] = myLDA(DD,G);
[lambda] = mynLDA(DD,G,1,3);

h = [];
p = [];
for l = 1:size(lambda,2)
    wv = [uC(1:9,:);uC(widx+18,:)]*lambda(:,l);
    mv = [uC(10:18,:);uC(midx+18,:)]*lambda(:,l);
    [h(l) p(l)] = ttest2(wv,mv);
    figure;
    
ksdensity(wv);
hold all
ksdensity(mv);
end

%close all
%{
splitVec1 = PCA_BKPROJ(lambda',uE1,uU1);
splitVec2 = PCA_BKPROJ(lambda',uE1,uU1+j1);
%}
SIG1 = bsxfun(@plus,uS1,j1);
SIG2 = bsxfun(@plus,uS2,j2);
%{
for e = 1:size(SIG,1)
    SIG(e,:) = cwt(-SIG(e,:),[50],'gaus2');
end
%}
%SIG = imfilter(SIG,fspecial('average',[1 150]),'replicate');
%SIG = gS;
%SIG = nnS;
%SIG = D;

mU_root = mean(SIG2(midx,:));
mS_root = std(SIG2(midx,:),1,1)*numel(midx)^-.5;

mU_cot = mean(SIG1(10:18,:));
mS_cot = std(SIG1(10:18,:),1,1)*numel(midx)^-.5;

wU_root = mean(SIG2(widx,:));
wS_root = std(SIG2(widx,:),1,1)*numel(widx)^-.5;


wU_cot = mean(SIG1(1:9,:));
wS_cot = std(SIG1(1:9,:),1,1)*numel(widx)^-.5;

%close all
figure;
errorbar(mU_root,mS_root,'r');
hold on
errorbar(mU_cot,mS_cot,'m');


errorbar(wU_root,wS_root,'b');
errorbar(wU_cot,wS_cot,'c');
hold on
%{
plot(.2*bsigL,'r');
plot(.2*bsigG1,'b')
plot(.2*bsigG2,'g');
%}
%{
CL = {'g' 'k' 'y' 'g--' 'k--'}
for l = 1:size(splitVec,1)
    if l ~= 2
        plot(2*(splitVec(l,:)-U),CL{l});
    else
        plot(-2*(splitVec(l,:)-U),CL{l});
    end
end
%}
LEG{1} = 'mut*root';
LEG{2} = 'mut*cot';
LEG{3} = 'wt*root';
LEG{4} = 'wt*cot';
legend(LEG);
%% take apart

%close all
bsigL = ((splitVec(1,:)-uU)) < 0;
bsigG = (splitVec(1,:)-uU) > 0;
bsigL = bwlarge(bsigL);
bsigG1 = bwlarge(bsigG);
bsigG2 = bwlarge(bsigG - bsigG1);

plot((splitVec(1,:)-uU));
hold on
plot(.2*bsigL,'r');
plot(.2*bsigG1,'b')
plot(.2*bsigG2,'g');
fidx1 = round(sum((1:numel(bsigL)).*bsigL)/sum(bsigL));
fidx2 = round(sum((1:numel(bsigG1)).*bsigG1)/sum(bsigG1));
fidx3 = round(sum((1:numel(bsigG2)).*bsigG2)/sum(bsigG2));
Z = zeros(size(bsigL));
%{
Z(fidx1) = -1;
Z(fidx3) = 1;
Z(fidx2) = 1;
%}
Z(112) =-1;
Z(64) = 1;
Z(32) = 1;
Z(end) = -1;
Z(191) = 1;
%Z = imfilter(Z,fspecial('gaussian',[1 50],10));
figure;

%Z(1:45) = splitVec(1,1:45) - uU(1:45)
%Z(273:end) = splitVec(1,273:end)-uU(273:end);
plot(Z)
hold on;
plot(splitVec(1,:)-uU);
%{
Z = splitVec(1,:)-uU;
Z(1:45) = 0;
Z(273:end) = 0;
%}
ssig = (splitVec(1,:)-uU);
ssig = Z;
figure;
plot(ssig);
%def = bsxfun(@minus,uS,uU)*ssig';
def = uS*ssig';
[th tp] = ttest2(def(widx+18),def(midx+18))
%def = uC*lambda(:,1);
[f1,x1] = ksdensity(def(widx +18));
%[f1,x1] = ksdensity(def(1:9));
figure;
hold on
plot(x1,f1,'k')
[f2,x2] = ksdensity(def(midx+18));
%[f2,x2] = ksdensity(def(10:18));
plot(x2,f2,'r')
%% lambda sweep

figure
shoot = 1;



if shoot
    tmpC = mean(uC1,1);
else
    tmpC = mean(uC2,1);
end

ti = 1;



%LAM = uC([1:9 widx 10:18 midx],:)*lambda(:,1);
%LAM = uC([1:18],:)*lambda(:,ti);
if shoot
    LAM = uC1*lambda(:,ti);
else
    LAM = uC2*lambda(:,ti);
end
%l = linspace(min(LAM),max(LAM),5);


if shoot
    LAM1 = LAM(1:9,:);
    LAM2 = LAM(10:end,:);
else
    LAM1 = LAM(widx,:);
    LAM2 = LAM(midx,:);
end
ul1 = mean(LAM1);
ul2 = mean(LAM2);
sl1 = .5*std(LAM1);
sl2 = .5*std(LAM2);
l = linspace(ul1-sl1,ul1+sl1,5);
for e = 1:numel(l)
    
    
    %sc = tmpC - lambda(:,ti)'*(tmpC*lambda(:,ti));
    sc = lambda(:,ti)'*(tmpC*lambda(:,ti));
    sc = sc + lambda(:,ti)'*l(e);
    if shoot
        % for cots
        sc = PCA_BKPROJ(sc,uE1,uU1+j1);
    else
        sc = PCA_BKPROJ(sc,uE2,uU2+j2);
    end
    % for roots
    
    plot(sc,'k');
    hold on
end


l = linspace(min(LAM2),max(LAM2),4);
l = linspace(ul2-sl2,ul2+sl2,5);
for e = 1:numel(l)
    sc = tmpC - lambda(:,ti)'*(tmpC*lambda(:,ti));
    sc = sc + lambda(:,ti)'*l(e);
    if shoot
        % for cots
        sc = PCA_BKPROJ(sc,uE1,uU1+j1);
    else
        sc = PCA_BKPROJ(sc,uE2,uU2+j2);
    end
    plot(sc,'r');
    hold on
end
if shoot
    sc = PCA_BKPROJ(tmpC,uE1,uU1+j1);
else
    sc = bsxfun(@plus,PCA_BKPROJ(tmpC,uE2,uU2),j2);
end
    plot(sc,'b');
    %plot(mean(nD,1),'y')
%% list names
for e = 1:numel(data)
    data(e).n
end
%% parse names - ver 1
prop = [];
for e = 1:numel(data)
    %[pth,nm,ext] = fileparts(data(e).n);
    nm = lower(data(e).n);
    prop(e,1) = ~(~isempty(strfind(lower(nm),'aca')) | ~isempty(strfind(lower(nm),'230')) | ~isempty(strfind(lower(nm),'231')) );
    prop(e,2) = ~isempty(strfind(nm,'chitin'));
    %prop(e,3) = isempty(strfind(nm,'chitin'));
    prop(e,3) = isempty(strfind(nm,'b222'));
end
UQ = unique(prop,'rows');
for e = 1:size(UQ,1)
    idx = all(bsxfun(@eq,UQ(e,:),prop),2);
    grp(e).idx = idx;
    grp(e).cnt = numel(idx);
    grp(e).grpLabel = UQ(e,:);
end
cmp(1,:) = [1 2];
cmp(2,:) = [1 3];
cmp(3,:) = [1 4];
cmp(4,:) = [1 5];
cmp(5,:) = [1 6];
cmp(6,:) = [3 6];

LAB{1,1} = 'WT';
LAB{2,1} = 'MUT';
LAB{1,2} = 'FLAG22';
LAB{2,2} = 'CHITIN';
LAB{1,3} = 'B222';
LAB{2,3} = 'NOTB222';
toPlot = S;
%toPlot = gS;
%toPlot = kS;;

values = sort(toPlot(:));
mx = mean(values(end-50:end));
close all
for e = 1:size(cmp,1)
    s = [];
    u = [];
    for g = 1:size(cmp,2)
        fidx = grp(cmp(e,g)).idx;
        u(g,:) = mean(toPlot(fidx,:),1);
        s(g,:) = std(toPlot(fidx,:),1,1);
        s(g,:) = s(g,:) * sum(fidx)^-.5;
        tmpL = [];
        for lab = 1:size(prop,2)
           tmpL = [tmpL  LAB{grp(cmp(e,g)).grpLabel(lab)+1,lab} '-'];
        end
        tmpL(end) = [];
        LABEL{g} = tmpL;
        errorbar(u(g,:),s(g,:));
        axis([0 320 -.1 mx])
        hold on
        title(LABEL);
        
        
    end
    legend(LABEL)
    %DELTA(e,:) = u(2,:) - u(1,:);
    hold off
    waitforbuttonpress
end

%{
close all
for e = 1:size(cmp,1)
    s = [];
    u = [];
    MM = [];
    for g = 1:size(cmp,2)
        fidx = grp(cmp(e,g)).idx;
        u = mean(toPlot(:,:,fidx),3);
        tmpL = [];
        for lab = 1:size(prop,2)
           tmpL = [tmpL  LAB{grp(cmp(e,g)).grpLabel(lab)+1,lab} '-'];
        end
        tmpL(end) = [];
        LABEL{g} = tmpL;
       
        MM = [MM u];
        %axis([0 320 -.1 mx])
        %hold on
        title(LABEL);
        
        
    end
    mesh(MM)
    legend(LABEL)
    %DELTA(e,:) = u(2,:) - u(1,:);
    hold off
    waitforbuttonpress
end
%}


%%
for e = 1:numel(data)
    sz(e) = numel(data(e).d);
end
%%
for e = 1:numel(data)
    data(e).d = data(e).d(1:min(sz));
end

%%
close all
for e = 1:size(uS,1)
    plot(uuD(e,:),'k');
    hold on
    plot(uS(e,:),'r');
    drawnow
    pause(.5);
    %waitforbuttonpress
    hold off
end
%% sweep data
close all
uC = mean(C,1);
for d = 1:size(C,2)
    l = linspace(min(C(:,d)),max(C(:,d)),7);
    for e = 1:numel(l)
        tmp = uC;
        tmp(d) = l(e);
        M = PCA_BKPROJ(tmp,E,U);
        plot(M)
        drawnow
        hold on
    end
    waitforbuttonpress
    close all
end