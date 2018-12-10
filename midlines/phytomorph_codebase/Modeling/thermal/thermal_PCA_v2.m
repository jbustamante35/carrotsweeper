gd = readtext('/home/nate/Downloads/myGDtemp.csv');
fprintf(['done reading dna markers']);
%%
md = cell2mat(gd(2:end,3:end));
%% my mappping
[mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL(md,size(md,2))
%% read in data file
%D = readtext('/home/nate/widef.csv');
D2 = readtext('/home/nate/Downloads/widef.csv');
fprintf(['done reading phenotypic data']);
%%
D1 = readtext('/home/nate/Downloadclos/widef_licor_bad_data_purged.csv');
%% select out the data from data sheet - NON-licor
X2 = D2(2:end,4:end);
%X = DL(2:end,(6:end-100));
% fill in empty with zeros
for e = 1:numel(X2)
    if isempty(X2{e})
        X2{e} = 0;
    end
end
% convert cell to mat
X2 = cell2mat(X2);
% transpose
X2 = X2';
% select and remove data
ridx = any(X2 == 0);
kidx = ~any(X2 == 0);
X2(:,find(ridx)) = [];
%% select out the data from data sheet - licor
X = D(2:end,5:end);
%X = DL(2:end,(6:end-100));
% fill in empty with zeros
for e = 1:numel(X)
    if isempty(X{e})
        X{e} = 0;
    end
end
% convert cell to mat
X = cell2mat(X);
% transpose
X = X';
% select and remove data
ridx = any(X == 0);
kidx = ~any(X == 0);
X(:,find(ridx)) = [];
%% find genotype
%{
clear G
for e = 2:size(D,1)
    fidx = strfind(D{e,1},' ');
    G(e-1) = str2num(D{e,1}(1:(fidx(1)-1)));
end
%}
G2 = D2(2:end,2);
G2(find(ridx)) = [];
%% build model for the NON-licor
[preVARE_X,preVARE_Y,preVARE_LX,preVARE_LY,preX_model,norX_model,preY_model,X2Y] = pre_post_model(X2,[240 241],3,[3 3]);

% percent predicted and not-predicted in Y
preVARE_Y
%% build model for the licor
[preVARE_X,preVARE_Y,preVARE_LX,preVARE_LY,preX_model,norX_model,preY_model,X2Y] = pre_post_model(X,[121 122],3,[3 3]);

% percent predicted and not-predicted in Y
preVARE_Y
%% build model
%[preVARE_X,preVARE_Y,preVARE_LX,preVARE_LY,preX_model,norX_model,preY_model] = pre_post_model(X,120,3,[3 3]);
%% model sweeper
close all
[paraPRE_X,paraNPRE_X,paraPRE_Y,preData,postData] = extractTraits(X,240,preX_model,norX_model,preY_model);
oPath = '/home/nate/Downloads/stomataData/';
mkdir(oPath);
modelSweeper(preX_model,'predictable-X',oPath);
modelSweeper(norX_model,'NONpredictable-X',oPath);
modelSweeper(preY_model,'NONpredictable-Y',oPath);
[npostS npostC npostU npostE npostL npostERR npostLAM] = PCA_FIT_FULL_T(postData,3);
%%
toSPLIT = paraPRE_Y;
%toSPLIT = paraPRE_X;
tmpModel = preY_model;
%tmpModel = preX_model;


% TURN ON SORTING!!!!
nG2 = [];
for e = 1:numel(G2)
    nG2{e} = (find(strcmp(G2{e},unique(G2))));
end
[~,sidx] = sort(cell2mat(nG2));
[~,isidx] = sort(sidx);


tmpC = nG2;
nG2 = nG2(sidx);
toSPLIT = toSPLIT(:,sidx);



mLDA = LDA(toSPLIT', nG2);
mLDA.Compute();


%TF = mLDA.Transform(paraPRE_Y', 1);
%dimension of a samples is < (mLDA.NumberOfClasses-1) so following line cannot be executed:
%transformedSamples = mLDA.Transform(meas, mLDA.NumberOfClasses - 1);

%transformedTrainSamples = mLDA.Transform(trainSamples, 1);

%{
UQ = unique(G2);
C = zeros(size(grandU,1));
%{
for u = 1:numel(UQ)
    fidx = find(strcmp(G2,UQ{u}));
    subD2(:,u) = mean(toSPLIT(:,fidx),2);
end
%}
grandU = mean(subD2,2);


for u = 1:numel(UQ)
    %fidx = find(G==UQ(u));
    fidx = find(strcmp(G2,UQ{u}));
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
G2 = G2(isidx);


%{



    k = gd(2:end,2);
    for e = 1:numel(k)
        k{e} = k{e}(2:end-1);
    end


    % align 
    UQ = unique(G2);
    GM = [];
    for u = 1:numel(k)
        fidx = find(strcmp(G2,k{u}));
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
    md = cell2mat(gd(2:end,3:end));
   
    uX = mean(md,1);
    uY = mean(GM,1);
    centerX = bsxfun(@minus,md,uX);
    centerY = bsxfun(@minus,GM,uY);
    options = statset('UseParallel',true,'Display','iter');
    [B FitInfo] = lassoglm(centerX(:,1:5250),centerY(:,1),'normal','CV',10,'Alpha',.4,'Options',options);
    


    [mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL(md,size(md,2));
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
     fidx = find(strcmp(G2,UQ{u}));
     tmpSub = toF(:,fidx);
     meanGG(:,u) = mean(tmpSub,2);
     varGG(:,u) = std(tmpSub,1,2).^2;
end
(std(meanGG,1,2).^2).*mean(varGG,2).^-1


toF = newH;
for u = 1:numel(UQ)
     %fidx = find(G==UQ(u));
     fidx = find(strcmp(G2,UQ{u}));
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
newD1{1,1} = D2{1,1};
newD1{1,2} = D2{1,2};
newD1{1,3} = D2{1,3};
for e = 1:size(paraPRE_X,1)
    newD1{1,3+e} = ['trait' num2str(e)];
end
for e = 1:(size(D2,1)-1)
    for k = 1:3
        newD1{1+e,k} = D2{1+e,k};
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
D2 = readtext('/home/nate/widef_licor.csv');

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


