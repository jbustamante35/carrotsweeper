%%%%%%%%% GENOTYPE %%%%%%%%%
%% read dna markers
rawGenotypeData = readtext('/home/nate/Downloads/myGDtemp_maf_05.csv');
fprintf(['done reading dna markers']);
%% get marker data
G = cell2mat(rawGenotypeData(2:end,3:end));
fprintf(['done section marker']);
%% get the chromosome locations
Gname = rawGenotypeData(1,3:end);
for e = 1:numel(Gname)
    Gname{e} = Gname{e}(2:end-1);
    fidx = strfind(Gname{e},'_');
    POS(e) = str2num(Gname{e}((fidx(end)+1):end));
end
fprintf(['done section chromosome position']);
%% get genotype names
GN = rawGenotypeData(2:end,2);
for e = 1:numel(GN)
    GN{e} = GN{e}(2:end-1);
end
%%%%%%%%% PHENOTYPE %%%%%%%%%
%% read in data file
rawPhenotypeData = readtext('/home/nate/Downloads/widef_covariates.csv');
fprintf(['done reading phenotypic data']);
%% select out the data from data sheet 
P = rawPhenotypeData(2:end,9:end)';
COP = cell2mat(rawPhenotypeData(2:end,3:7));
% fill in empty with zeros
for e = 1:numel(P)
    if isempty(P{e})
        P{e} = 0;
    end
end
% convert cell to mat
P = cell2mat(P);
% transpose
P = P';
% select and remove data
ridx = any(P == 0);
kidx = ~any(P == 0);
P(:,find(ridx)) = [];
%% find genotype
PN = rawPhenotypeData(2:end,2);
PN(find(ridx)) = [];
%% intersect between GN and PN
iN = intersect(PN,GN);
rm = [];
for e = 1:numel(PN)
    if ~isempty(find(strcmp(PN{e},GN)))
        rm(e) = 0;
    else
        rm(e) = 1;
    end
end
PN(find(rm==1)) = [];
rm = [];
for e = 1:numel(GN)
    if ~isempty(find(strcmp(GN{e},PN)))
        rm(e) = 0;
    else
        rm(e) = 1;
    end
end
GN(find(rm==1)) = [];
%% correct for pop structure
[genoS genoC genoU genoE genoL genoERR genoLAM] = PCA_FIT_FULL(G,4);
%%
P = mean(P(:,(end-20):end),2);
%% COP 
copC1 = [ones(size(COP,1),1) COP];
%copC1 = COP;
%{
netCOP = newgrnn(copC1',P');



copBeta = copC1\P;
%}
net = feedforwardnet([5 5]);
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(P,4);
netCOP = train(net,copC1',pC','UseParallel','yes');
copP = sim(netCOP,copC1')';
copP = PCA_BKPROJ(copP,pE,pU);


for r = 1:100
    plot(P(r,:),'k');
    hold all
    plot(copP(r,:),'r');
    axis([0 size(P,2) -.1 1])
    waitforbuttonpress
    hold off
end
%%
P = P - copP;
%% build model based on pc
PU_level1 = [];
P2 = [];
km = 0;
for e = 1:numel(GN)
    fidx = find(strcmp(GN{e},PN));
    PU_level1(e,:) = mean(P(fidx,:),1);
end
genoC1 = [ones(size(genoC,1),1) genoC];
%{

%genoC1 = [genoC];
genoBeta = genoC1\PU_level1;


net = newgrnn(genoC1',PU_level1');
%net = feedforwardnet([4 2]);
%net = train(net,genoC1',PU_level1');
preP = sim(net,genoC1')';

if km
    kidx = kmeans(preP,2);
    net2 = newgrnn(kidx',preP');
    preP2 = sim(net2,kidx')';
end
%}
[kidx] = kmeans(genoC, 20);
UQ = unique(kidx);
for u = 1:numel(UQ)
    Gu(u,:) = mean(genoC(kidx==u,:),1);
    Pu(u,:) = mean(P(kidx==u,:),1);
    plot3(genoC(kidx==u,1),genoC(kidx==u,2),genoC(kidx==u,3),'.')
    hold all
end
net = newgrnn(Gu',Pu',.0001);
preP = sim(net,genoC')';


%{
net = feedforwardnet([5 5]);
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(PU_level1,4);
netCOP = train(net,genoC1',pC','UseParallel','yes');
copP = sim(netCOP,genoC1')';
preP = PCA_BKPROJ(copP,pE,pU);
%}

for r = 1:10
    plot(P(r,:),'k');
    hold all
    plot(preP(r,:),'r');
    waitforbuttonpress
    hold off
end
%preP = genoC1*genoBeta;
mean(sum((preP - PU_level1).^2,2))
for e = 1:numel(PN)
    fidx = find(strcmp(PN{e},GN));
    P2(e,:) = bsxfun(@minus,P(fidx,:),preP(fidx,:));
end
%% build model for the NON-licor
[preVARE_X,preVARE_Y,preVARE_LX,preVARE_LY,preX_model,norX_model,preY_model,X2Y] = pre_post_model(P',[240 241],3,[3 3]);

% percent predicted and not-predicted in Y
preVARE_Y
%% model sweeper
close all
[paraPRE_X,paraNPRE_X,paraPRE_Y,preData,postData] = extractTraits(P',240,preX_model,norX_model,preY_model);
oPath = '/home/nate/Downloads/stomataData/';
mkdir(oPath);
modelSweeper(preX_model,'predictable-X',oPath);
modelSweeper(norX_model,'NONpredictable-X',oPath);
modelSweeper(preY_model,'NONpredictable-Y',oPath);
[npostS npostC npostU npostE npostL npostERR npostLAM] = PCA_FIT_FULL_T(postData,3);
%{

%%
toSPLIT = paraPRE_Y;
%toSPLIT = paraPRE_X;
tmpModel = preY_model;
%tmpModel = preX_model;

% sort the data
% TURN ON SORTING!!!!
nG2 = [];
for e = 1:numel(PN)
    nG2{e} = (find(strcmp(PN{e},unique(PN))));
end
[~,sidx] = sort(cell2mat(nG2));
[~,isidx] = sort(sidx);
% sort with sidx
nG2 = nG2(sidx);
toSPLIT = toSPLIT(:,sidx);
%}
%% match for FUL
toSPLIT = paraPRE_Y';
colGN = [];
for e = 1:numel(PN)
    fidx = (find(strcmp(PN{e},GN)));
    fidx
    if ~isempty(fidx)
        colGN(e) = fidx;
    else
        colGN(e) = 0;
    end
end
rm = colGN == 0;
colGN(rm) = [];
toSPLIT(rm,:) = [];
fD = postData';
fD(rm,:) = [];
eD = preData';
eD(rm,:) = [];
%% make geno type means
for e = 1:numel(GN)
    fidx = find(strcmp(GN{e},PN));
    uP(e,:) = mean(toSPLIT(fidx,:),1);
end
%%
pG = G';
pG = reshape(pG,[size(pG,1) 1 1 size(pG,2)]);
N = 15000;
layers = [ ...
    imageInputLayer([N 1 1])
    convolution2dLayer([23 1],3)
    reluLayer
    maxPooling2dLayer([10 1],'Stride',1)
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','InitialLearnRate',0.001,'MaxEpochs',300,'ExecutionEnvironment','parallel');
C = [];
for tr = 1:10
    [Train, Test] = crossvalind('HoldOut', 661, .2);
    net = trainNetwork(pG(1:N,:,:,Train),uP(Train,1),layers,options);
    preC = predict(net,pG(1:N,:,:,Test));
    C(tr) = corr(preC,uP(Test,1));
    C
end

%% MY LDA SCAN
NT = size(G,2);
lambda = zeros(3,NT);
for c = 1:NT
    tic
        grp1 = G(:,c) == 0;% | G(:,c) == 1;
        grp2 = ~grp1;


        grp1 = find(grp1);
        grp2 = find(grp2);

        LAB = zeros(size(toSPLIT,1),1);
        for i = 1:numel(grp2)
            LAB(grp2(i) == colGN) = 1;
        end
        NUM(c,:) = [numel(grp1) numel(grp2)];
        if ~isempty(grp1) & ~isempty(grp2)
            lambda(:,c) = myLDA(toSPLIT,LAB);
            tmpC = toSPLIT*lambda(:,c);
            %[h(c) p(c)] = ttest2(tmpC(LAB==0),tmpC(LAB==1));
            p(c) = anova1(tmpC,LAB,'off');
        else
            lambda(:,c) = 0;
            h(c) = 0;
            p(c) = inf;
        end
    toc
    c
    size(G,2)
end
%% clip out chro
close all
for e = 1:numel(Gname)
    if ~isempty(strfind(Gname{e},'pstI__1'))
    %if ~isempty(strfind(Gname{e},'"pstI__3_66050'))
    %if ~isempty(strfind(Gname{e},'52156808'))
    %if ~isempty(strfind(Gname{e},'52017857'))
        PUL(e) = 1;
    else
        PUL(e) = 0;
    end
end
disp = 1;
% fun fun fun final
TOT = [];
TOTM = [];
TOTP = [];
TOTQ = [];
TOTF = [];
for pp = 1:size(pTEST,2)
    toPLOT = pTEST(:,pp);
    
    
    
    [sp,spi] = sort(toPLOT,'ascend');
    [~,ispi] = sort(spi);
    CON = sum((1:numel(p)).^-1);
    %qSEL = .01;
    %CON = 1;
    numel(p)*qSEL
    wow = (1:numel(p))*(numel(p)^-1)*(qSEL).*CON.^-1;
    kidx = (sp < wow');
    WW = sp < .05;
    plot(kidx.*WW)
    sum((kidx.*WW))
    k = find(kidx);
    THRES = sp(k(end));
    
    pp
    if disp
        find(PUL);
        close all
        plot(POS(find(PUL)),-log10(toPLOT(find(PUL))));
        tmpF = find(PUL);
        [pv,tmpI] = min(toPLOT(find(PUL)));
        hold on
        POS(tmpF(tmpI))
        plot(POS(tmpF(tmpI)),-log10(toPLOT(tmpF(tmpI))),'r*')
        %plot(POS(midx(ssidx(1))),-log10(toPLOT(midx(ssidx(1)))),'g*')
        plot(POS(find(PUL)),-log10(.05/numel(toPLOT))*ones(1,numel(find(PUL))),'r');
    end
    
    
    [FDR,Q,P10] = mafdr(toPLOT);
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(toPLOT,.05);
    %[pFDR] = mafdr(toPLOT,'BHFDR',true);
    TOT(pp) = sum(FDR < .0000001);
    TOTM(pp) = sum(FDR < Q);
    TOTQ(pp) = sum(Q < .00001);
    TOTP(pp) = sum(toPLOT < THRES);
    TOTF(pp) = sum(pFDR < .05);
    figure;
    plot(POS(find(PUL)),(-log10(toPLOT(find(PUL)))));
    hold on
    plot(POS(find(PUL)),-log10((THRES))*ones(1,numel(find(PUL))),'g');
    %waitforbuttonpress
    
    
    figure;
    plot(POS(find(PUL)),(-log10(FDR(find(PUL)))));
    hold on
    plot(POS(find(PUL)),-log10((.0001))*ones(1,numel(find(PUL))),'g');
   % waitforbuttonpress
end
numel(toPLOT)*.0000001
%%

%%
close all
qSEL = .1;
for qSEL = linspace(.05,.0001,10)
    p = pTEST(:,1);
    [sp,spi] = sort(p,'ascend');
    [~,ispi] = sort(spi);
    CON = sum((1:numel(p)).^-1);
    %qSEL = .01;
    %CON = 1;
    numel(p)*qSEL
    wow = (1:numel(p))*(numel(p)^-1)*(qSEL).*CON.^-1;
    kidx = (sp < wow');
    WW = sp < .05;
    plot(kidx.*WW)
    sum((kidx.*WW))
    k = find(kidx);
    THRES = sp(k(end));
    hold on
    waitforbuttonpress
end
%%  apply 
%plot(proj);
aVEC = 1:10;
proj = toSPLIT*lambda(:,spi(1));
proj = mean(postData(end-20:end,:)',2);
%proj = P2;

%for v = 1:numel(aVEC)
    %proj = toSPLIT*lambda(:,spi(1:10));
    hTEST = zeros(NT,size(proj,2));
    pTEST = zeros(NT,size(proj,2));
    for pr = 1:100
        pPROJ(:,pr) = proj(randperm(size(proj,1)),1);
    end
    for c = 1:NT
        tic
            grp1 = G(:,c) == 0;% | G(:,c) == 1;
            grp2 = ~grp1;


            grp1 = find(grp1);
            grp2 = find(grp2);

            LAB = zeros(size(toSPLIT,1),1);
            for i = 1:numel(grp2)
                LAB(grp2(i) == colGN) = 1;
            end
    
            
            
            if ~isempty(grp1) & ~isempty(grp2)
                for v = 1:size(proj,2)
                    %[hTEST(c,v) pTEST(c,v)] = ttest2(proj(LAB==0,v),proj(LAB==1,v));
                    pTEST(c,v) = anova1(proj(:,v),LAB,'off');
                end
                %{
                for v = 1:size(proj,2)
                    for pr = 1:10
                        %[J,tmpJ(pr)] = ttest2(pPROJ(LAB==0,pr),pPROJ(LAB==1,pr));
                        tmpJ(pr) = anova1(pPROJ(:,pr),LAB,'off');
                    end
                    pRTEST(c,v) = mean(tmpJ);
                end
                %}
                
            else
                hTEST(c,:) = 0;
                pTEST(c,:) = inf;
            end
        toc
        c
        size(G,2)
    end
%end
%%

%%
NUM(spi(1:40),:)

[mv,midx] = min(p);
s = preY_model.sim(lambda(:,midx)) - preY_model.U;
figure;
plot(s)
plot(X2Y(mean(eD(LAB==0,:))'),'k')
hold on
plot(X2Y(mean(eD(LAB==1,:))'),'r')
plot(X2Y(mean(eD(LAB==2,:))'),'b')
%%
[~,MMidx] = min(pTEST,[],1);
close all
%c = MMidx(2);
[~,c] = min(p);
Gname(c)
[~,phune] = sort(p);
grp1 = G(:,c) == 0;% | G(:,c) == 1;
grp2 = G(:,c) == 1;
grp3 = G(:,c) == 2;

grp1 = find(grp1);
grp2 = find(grp2);
grp3 = find(grp3);

LAB = zeros(size(toSPLIT,1),1);

for i = 1:numel(grp2)
        LAB(grp2(i) == colGN) = 1;
end
    
for i = 1:numel(grp3)
        LAB(grp3(i) == colGN) = 2;
end
%{
close all
midx = find(NUM(:,1) < 350 & NUM(:,1) > 310 & p' < .01);
midx = intersect(midx,find(PUL));
[spp ssidx] = sort(p(midx));
midx(ssidx(1))
%}
c1 = preY_model.func(mean(eD(LAB==0,:))',mean(fD(LAB==0,:))');
s1 = preY_model.sim(c1);

c2 = preY_model.func(mean(eD(LAB==1,:))',mean(fD(LAB==1,:))');
s2 = preY_model.sim(c2);

c3 = preY_model.func(mean(eD(LAB==02,:))',mean(fD(LAB==2,:))');
s3 = preY_model.sim(c3);
figure;
plot(s1,'r');hold all;plot(s2,'g');plot(s3,'b');

figure;
plot(mean(fD(LAB==0,:)));hold all;plot(mean(fD(LAB==1,:)));plot(mean(fD(LAB==2,:)));

figure;
errorbar(mean(fD(LAB==0,:)),std(fD(LAB==0,:),1,1)*sum(LAB==0)^-.5,'r');hold all;
errorbar(mean(fD(LAB==1,:)),std(fD(LAB==1,:),1,1)*sum(LAB==1)^-.5,'g');
errorbar(mean(fD(LAB==2,:)),std(fD(LAB==2,:),1,1)*sum(LAB==2)^-.5,'b');

%%
plot(s);
s = preY_model.sim(lambda(:,e));
figure;
plot(X2Y(mean(eD(LAB==0,:))'),'k')
hold on
plot(X2Y(mean(eD(LAB==1,:))'),'r')


%%
close all
for e = 1:size(lambda,2)
    s = preY_model.sim(lambda(:,e)) - preY_model.U;
    fs = imfilter(s,fspecial('average',21),'replicate');
    zidx = find(fs(1:end-1) > 0 & fs(2:end) < 0 | fs(1:end-1) < 0 & fs(2:end) > 0);
    nZ(e) = numel(zidx);
    Z = zeros(size(s));
    Z(zidx) = 1;
    %{
    plot(s)
    hold on
    plot(Z);
    axis([0 350 -3*10^-5 3*10^-5])
    hold off
    drawnow
    %hold on
    %}
    e
end

%% group via zero crossing
close all
histogram(nZ)
i2 = find(nZ==2);
s = bsxfun(@minus,preY_model.sim(lambda(:,i2)),preY_model.U);
for e = 1:size(s,2)
    s(:,e) = s(:,e)/norm(s(:,e));
end
seedN = 100;
subI = [];
SN = [];
for e = 1:size(s,2)
    if abs(s(:,e)'*s(:,seedN)) > .8
        subI = [subI e];
        if s(:,e)'*s(:,seedN) < 0
            SN = [SN -1];
        else
            SN = [SN 1];
        end
    end
end
plot(s(:,subI));
sd = setdiff(1:numel(i2),subI);
%%
close all
ss = bsxfun(@times,s(:,subI),SN);
tmpLambda = bsxfun(@times,lambda(:,i2(subI)),SN);
for e = 1:size(tmpLambda,2)
    tmpLambda(:,e) = tmpLambda(:,e) / norm(tmpLambda(:,e));
end
[lS lC lU lE lL lERR lLAM] = PCA_FIT_FULL_T(ss,2);
[lS lC lU lE lL lERR lLAM] = PCA_FIT_FULL_T(tmpLambda,2);
%sweep

for g = 1:size(lC,1)
    tmpUU = mean(lC,2);
    LL = linspace(min(lC(g,:)),max(lC(g,:)),5);
    for l = 1:numel(LL)
        tmpUU(g) = LL(l);
        M = PCA_BKPROJ_T(tmpUU,lE,lU);
        plot(M)
        hold on
    end
    waitforbuttonpress
    hold off
end
myL = [];
myL(:,1) = PCA_BKPROJ_T([min(lC(1,:));0],lE,lU);
myL(:,2) = PCA_BKPROJ_T([max(lC(1,:));0],lE,lU);
myL(:,3) = PCA_BKPROJ_T([0;0],lE,lU);
myL(:,4) = PCA_BKPROJ_T([0;min(lC(2,:))],lE,lU);
myL(:,5) = PCA_BKPROJ_T([0;max(lC(2,:))],lE,lU);

myLS = bsxfun(@minus,preY_model.sim(myL),preY_model.U);

%myL = bsxfun(@plus,myL,preY_model.U);
%%
mLDA = LDA(toSPLIT', nG2);
mLDA.Compute();


%TF = mLDA.Transform(paraPRE_Y', 1);
%dimension of a samples is < (mLDA.NumberOfClasses-1) so following line cannot be executed:
%transformedSamples = mLDA.Transform(meas, mLDA.NumberOfClasses - 1);

%transformedTrainSamples = mLDA.Transform(trainSamples, 1);

%{
UQ = unique(PN);
C = zeros(size(grandU,1));
%{
for u = 1:numel(UQ)
    fidx = find(strcmp(PN,UQ{u}));
    subD2(:,u) = mean(toSPLIT(:,fidx),2);
end
%}
grandU = mean(subD2,2);


for u = 1:numel(UQ)
    %fidx = find(G==UQ(u));
    fidx = find(strcmp(PN,UQ{u}));
    subD = toSPLIT(:,fidx);
    tmpU = mean(subD,2);
    tmpU = (tmpU - grandU);
    C = C + tmpU*tmpU';
    toSPLIT(:,fidx) = bsxfun(@minus,toSPLIT(:,fidx),tmpU);
end
C = C * numel(UQ)^-1;
[V,DI] = eigs(inv(cov(toSPLIT'))*C);


newH = (toSPLIT'*V)';
%newH = newH(:,isidx);

%}

TF = mLDA.Transform(toSPLIT', 3);
V = mLDA.EigenVectors;
mLDA.EigenValues


EY = eye(3);
for e = 1:3
    top = EY(:,e)'*mLDA.BetweenScatter*EY(:,e);
    bot = EY(:,e)'*mLDA.WithinScatter*EY(:,e);
    R(e) = top / bot;
end
R'



newH = (toSPLIT'*V)';

tmpModel.E = tmpModel.E*V;
tmpModel.V = std(newH,1,2);
modelSweeper(tmpModel,'H-Model',oPath);

newH = newH(:,isidx);
PN = PN(isidx);


%{



    k = rawGenotypeData(2:end,2);
    for e = 1:numel(k)
        k{e} = k{e}(2:end-1);
    end


    % align 
    UQ = unique(PN);
    GM = [];
    for u = 1:numel(k)
        fidx = find(strcmp(PN,k{u}));
        subD = newH(:,fidx);
        tmpU = mean(subD,2)';
        GM(u,:) = tmpU;
    end

    fidx = 0;

    for u = 1:numel(UQ)
        fidx =  find(strcmp(UQ{u},k));
        if ~isempty(fidx)
            rs(u) = fidx;
        else
            rs(u) = 0;
        end
        fidx = 0;
    end
    rs(rs==0) = [];
    UQ(rs)
    G = cell2mat(rawGenotypeData(2:end,3:end));
   
    uX = mean(G,1);
    uY = mean(GM,1);
    centerX = bsxfun(@minus,G,uX);
    centerY = bsxfun(@minus,GM,uY);
    options = statset('UseParallel',true,'Display','iter');
    [B FitInfo] = lassoglm(centerX(:,1:5250),centerY(:,1),'normal','CV',10,'Alpha',.4,'Options',options);
    


    [mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL(G,size(G,2));
    lassoPlot(B,FitInfo,'plottype','CV');
    find(B(:,FitInfo.Index1SE))



%}
%{
mLDA = LDA(newH', nG2);
mLDA.Compute();
mLDA.EigenValues


EY = eye(3);
for e = 1:3
    top = EY(:,e)'*mLDA.BetweenScatter*EY(:,e);
    bot = EY(:,e)'*mLDA.WithinScatter*EY(:,e);
    R(e) = top / bot;
end
R'
%}
%{
toF = toSPLIT;
for u = 1:numel(UQ)
     %fidx = find(G==UQ(u));
     fidx = find(strcmp(PN,UQ{u}));
     tmpSub = toF(:,fidx);
     meanGG(:,u) = mean(tmpSub,2);
     varGG(:,u) = std(tmpSub,1,2).^2;
end
(std(meanGG,1,2).^2).*mean(varGG,2).^-1


toF = newH;
for u = 1:numel(UQ)
     %fidx = find(G==UQ(u));
     fidx = find(strcmp(PN,UQ{u}));
     tmpSub = toF(:,fidx);
     meanGG(:,u) = mean(tmpSub,2);
     varGG(:,u) = std(tmpSub,1,2).^2;
end
(std(meanGG,1,2).^2).*mean(varGG,2).^-1
%Hmodel = preY_model;
%}

%{

mLDAT1 = LDA(newH(1,:)', nG2);
mLDAT1.Compute();
mLDAT1.BetweenScatter*mLDAT1.WithinScatter^-1
mLDAT2 = LDA(newH(2,:)', nG2);
mLDAT2.Compute();
mLDAT2.BetweenScatter*mLDAT2.WithinScatter^-1
mLDAT3 = LDA(newH(3,:)', nG2);
mLDAT3.Compute();
mLDAT3.BetweenScatter*mLDAT3.WithinScatter^-1
%{
mLDATT = LDA(newH', nG2);
mLDATT.Compute();
mLDATT.EigenValues.^.5
%}
mLDAT1 = LDA(toSPLIT(1,:)', nG2);
mLDAT1.Compute();
mLDAT1.BetweenScatter/mLDAT1.WithinScatter^-1
mLDAT2 = LDA(toSPLIT(2,:)', nG2);
mLDAT2.Compute();
mLDAT2.BetweenScatter/mLDAT2.WithinScatter^-1
mLDAT3 = LDA(toSPLIT(3,:)', nG2);
mLDAT3.Compute();
mLDAT3.BetweenScatter/mLDAT3.WithinScatter^-1
%{
mLDATT = LDA(toSPLIT', nG2);
mLDATT.Compute();
mLDATT.EigenValues.^.5
%}

%}
%% make predicted sheet
close all
preY = X2Y(preData);
preDIF = preY_model.sim(paraPRE_Y);
numView = 30;
for e = 1:numView
    plot(preY(:,e),'r');
    hold on
    plot(postData(:,e),'k');
    plot(preDIF(:,e)+preY(:,e),'g');
    hold off
    axis([0 200 -.1 1]);
    drawnow
    %waitforbuttonpress
    pause(.1)
end
%% 
%% fill in the sheet information
clear newD1;
newD1{1,1} = rawPhenotypeData{1,1};
newD1{1,2} = rawPhenotypeData{1,2};
newD1{1,3} = rawPhenotypeData{1,3};
for e = 1:size(paraPRE_X,1)
    newD1{1,3+e} = ['trait' num2str(e)];
end
for e = 1:(size(rawPhenotypeData,1)-1)
    for k = 1:3
        newD1{1+e,k} = rawPhenotypeData{1+e,k};
    end
end
newD2 = newD1;
newD3 = newD1;
newD4 = newD1;
newD5 = newD1;
newD6 = newD1;
%% fill in newD1
fillWith = paraPRE_X;
fidx = find(kidx);
for e = 1:size(fillWith,2)
    for k = 1:size(fillWith,1)
        newD1{fidx(e)+1,3+k} = fillWith(k,e);
    end
end
%% fill in newD2
fillWith = paraNPRE_X;
fidx = find(kidx);
for e = 1:size(fillWith,2)
    for k = 1:size(fillWith,1)
        newD2{fidx(e)+1,3+k} = fillWith(k,e);
    end
end
%% fill in newD3
fillWith = paraPRE_Y;
fidx = find(kidx);
for e = 1:size(fillWith,2)
    for k = 1:size(fillWith,1)
        newD3{fidx(e)+1,3+k} = fillWith(k,e);
    end
end
%% fill in newD4
fillWith = newH;
fidx = find(kidx);
for e = 1:size(fillWith,2)
    for k = 1:size(fillWith,1)
        newD4{fidx(e)+1,3+k} = fillWith(k,e);
    end
end
%% newD5
fillWith = preY;
fidx = find(kidx);
for e = 1:size(fillWith,2)
    for k = 1:size(fillWith,1)
        newD5{fidx(e)+1,3+k} = fillWith(k,e);
    end
end
%% newD6
fillWith = preDIF;
fidx = find(kidx);
for e = 1:size(fillWith,2)
    for k = 1:size(fillWith,1)
        newD6{fidx(e)+1,3+k} = fillWith(k,e);
    end
end
%%
cell2csv('/mnt/snapper/nate/phenoType1.csv',newD1);
cell2csv('/mnt/snapper/nate/phenoType2.csv',newD2);
cell2csv('/mnt/snapper/nate/phenoType3.csv',newD3);
cell2csv('/mnt/snapper/nate/phenoTypeXH.csv',newD4);
cell2csv('/mnt/snapper/nate/phenoType5.csv',newD5);
cell2csv('/mnt/snapper/nate/phenoType6.csv',newD6);
%% make snp column data
F = readtext('/home/nate/Downloads/possible_QTL.csv');
GM = {};
toSCAN = 3;
for e = 2:size(F,1)
    GM{e-1,1} = F{e,2}(2:end-1);
    GM{e-1,2} = F{e,toSCAN};
end
for e = 1:size(GM,1)
    fidx = find(strcmp(G,GM{e,1}));
    NG(fidx) = GM{e,2};
end
% plot density
close all
wholeStack = [paraPRE_X;paraPRE_Y];
%wholeStack = npostC;
for e = 1:size(wholeStack,1)
    [f1 x1] = ksdensity(wholeStack(e,NG==0));
    [f2 x2] = ksdensity(wholeStack(e,NG==1));
    [f3 x3] = ksdensity(wholeStack(e,NG==2));
    
    [p,t,stats] = anova1(wholeStack(e,:),NG);
    multcompare(stats)
    plot(x1,f1,'k');
    hold on
    plot(x2,f2,'r');
    plot(x3,f3,'b');
    waitforbuttonpress
    hold off
    close all
end
%%
toSPLIT = wholeStack;
UQ = unique(NG);
grandU = mean(toSPLIT,2);
C = zeros(size(grandU,1));
for u = 1:numel(UQ)
    fidx = find(NG==UQ(u));
    %fidx = find(strcmp(G,UQ{u}));
    subD = toSPLIT(:,fidx);
    tmpU = mean(subD,2);
    tmpU = tmpU - grandU;
    C = C + tmpU*tmpU';
    toSPLIT(:,fidx) = bsxfun(@minus,toSPLIT(:,fidx),tmpU);
end
C = C * numel(UQ)^-1;
[V,DI] = eigs(inv(cov(toSPLIT'))*C);
splitPH = V'*toSPLIT;
%%
for e = 1:numel(NG)
    mNG{e} = NG(e);
end
mLDA = LDA(splitPH',mNG);
mLDA.Compute();
testF = real(mLDA.Transform(splitPH',3))';
%%
% plot density
close all
wholeStack = testF;
for e = 1:size(wholeStack,2)
    [f1 x1] = ksdensity(wholeStack(e,NG==0));
    [f2 x2] = ksdensity(wholeStack(e,NG==1));
    [f3 x3] = ksdensity(wholeStack(e,NG==2));
    plot(x1,f1,'k');
    hold on
    plot(x2,f2,'r');
    plot(x3,f3,'b');
    waitforbuttonpress
    hold off
end
%%
close all
toPlot = paraPRE_X;
%toPlot = paraNPRE_X;
toPlot = paraPRE_Y;
figure
ksdensity(toPlot(1,:))
figure;
ksdensity(toPlot(2,:))
figure
ksdensity(toPlot(3,:))
%%

%%

%% sweep models

%%
rawPhenotypeData = readtext('/home/nate/widef_licor.csv');

%%

%%
[S C U E L ERR LAM] = PCA_FIT_FULL_T(X,5);
%% view the data
close all
for e = 1:size(X,2)
    plot(X(:,e),'k');
    hold on
    plot(S(:,e),'r');
    hold off
    drawnow
    pause(.01);
end
%% sweep
close all
numToSweep = size(C,1);
numToSweep = 3;
for c = 1:numToSweep
    tmpU = mean(C,2);
    L = linspace(-std(C(c,:),1,2),std(C(c,:),1,2),7);
    for l = 1:numel(L)
        tmpC = tmpU;
        tmpC(c) = L(l);
        sweepC = PCA_BKPROJ_T(tmpC,E,U);
        plot(sweepC);
        hold all
    end
    waitforbuttonpress
    close all
end
%% ossiliations
front = X(1:400,:);
front = bsxfun(@minus,front,mean(front,1));
fT = fft(front,[],1);
[fS fC fU fE fL fERR fLAM] = PCA_FIT_FULL_T(abs(fT),3);
plot(abs(fT));
freq = 23/size(front,1);
T = freq^-1;
%%
[frS frC frU frE frL frERR frLAM] = PCA_FIT_FULL_T(front',10);
[frS frC frU frE frL frERR frLAM] = PCA_FIT_FULL_T(front',10);
%% correlations between freq and dynamcis
corr(fC',C');
%% fitting fun
back = X(401:end,:);
plot(back)
[bS bC bU bE bL bERR bLAM] = PCA_FIT_FULL_T(back,15);
[Xloadings,Yloadings,Xscores,Yscores, beta,pctVar,mse,stats,Weights] = plsregress(front',bC(:,:)',7);
[Xloadings,Yloadings,Xscores,Yscores, beta,pctVar,mse,stats,Weights] = plsregress(abs(fT)',bC(:,:)',7);
[Xloadings,Yloadings,Xscores,Yscores, beta,pctVar,mse,stats,Weights] = plsregress(front',back',5);
preD = ([ones(size(front',1),1) front']*beta)';
%plot(Yp(:,1),bC(1,:),'.');
%Yp = [ones(size(front',1),1) front']*beta;
%preD = PCA_BKPROJ_T(Yp',bE,bU);
%% non_predictive space
[wS wX wU wE wL wERR wLAM] = PCA_FIT_FULL_T(Weights,5);
% not in pre space
nF = bsxfun(@minus,front,mean(front,2));
tmpC = PCA_REPROJ_T(nF,wE,wU);
nF = PCA_BKPROJ_T(tmpC,wE,wU);
nF = bsxfun(@plus,nF,mean(front,2));

[nS nC nU nE nL nERR nLAM] = PCA_FIT_FULL_T(front - nF,3);


[aS aC aU aE aL aERR nLAM] = PCA_FIT_FULL_T(nF,5);

%% formal linear regress

%%%%%%%%%%%%%%%%%%%%%
close all
h1 = figure;
h2 = figure;
rX = front';
[jS rX jU jE jL jERR jLAM] = PCA_FIT_FULL(rX,5);
%rX = profileTOPX1000;
rY = back';
[jS rY jU jE jL jERR jLAM] = PCA_FIT_FULL(rY,5);
COR = [];
Yp = [];
for l = 1:10
    parfor e = 1:size(rX,1)
        idx = setdiff(1:size(rX,1),e);
        
        subX = rX(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        %size(R);
        subY = rY(idx,:);
        
        
        %{
        %subX = [Co(idx,:) LV(idx,1)];
        subX = [Co(idx,:)];
        %}
        net = feedforwardnet([l 5]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        tmp = sim(net,rX(e,:)');
        tmp = PCA_BKPROJ(tmp',jE,jU);
        
        Yp(l,e,:) = tmp;
        
        %{
        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
        %Yp = [ones(size(X,1),1) X(e,:)]*beta;
        %Yp(e) = [1 X(e,:)]*beta;
        %Yp(l,e) = [1 X(e,:) LV(e,1)]*beta;
        Yp(l,e,:) = [1 rX(e,:)]*beta;
        %}
        %{
        figure(h2)
        plot(rY(e),Yp(l,e),'.')
        drawnow
        hold on
        %}
        e
    end
    %R = corr(rY,squeeze(Yp(l,:,:)));
    R = corr(jS,squeeze(Yp(l,:,:)));
    COR(l,:) = diag(R);
    %COR(l) = corr(rY,Yp(l,:));
    %l
    figure(h1)
    plot(COR')
    drawnow
end
%% show predictions
close all
for e = 1:size(preD,2)
    plot(preD(:,e),'r')
    hold on
    plot(back(:,e),'k')
    hold off
    axis([0 600 0 .8]);
    drawnow
    pause(.01)
end
%{
%% normalized
dX = gradient(X')';
[dS dC dU E L ERR LAM] = PCA_FIT_FULL_T(dX,5);
%}
%%
P = back - preD;
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL_T(P,3);

close all
numToSweep = size(pC,1);
numToSweep = 3;

%{
sweepC = pC;
sweepE = pE;
sweepU = pU;
%}

sweepC = aC;
sweepE = aE;
sweepU = aU;


sweepC = nC;
sweepE = nE;
sweepU = nU;

for c = 1:numToSweep
    tmpU = mean(sweepC,2);
    L = linspace(-std(sweepC(c,:),1,2),std(sweepC(c,:),1,2),7);
    for l = 1:numel(L)
        tmpC = tmpU;
        tmpC(c) = L(l);
        sweepSHOW = PCA_BKPROJ_T(tmpC,sweepE,sweepU);
        
        
        plot(sweepSHOW);
        hold all
    end
    waitforbuttonpress
    close all
end

%%
for e = 1:size(P,2)
    sig = pS(:,e) + preD(:,e);
    plot(sig,'r');
    hold on
    plot(back(:,e),'k');
    hold off
    drawnow
    pause(.5)
end


