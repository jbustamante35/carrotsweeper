%%%%%%%%%%%%%%%%%%%%%%
% 1: run this to load data from cyverse
%%%%%%%%%%%%%%%%%%%%%
FilePath = '/mnt/spaldingdata/nate/mirror_images/carrotData/turnersarahd/return/wholeData_OLD/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
D = [];
D1 = [];
NM = {};
TOP = [];
BOT = [];
% [petCount newPetCount petWidth petW mean(petL) petNUM petDIA]
for e = 1:numel(FileList)
    [p,nm,ext] = fileparts(FileList{e});
    if ~isempty(strfind(nm,'topProfile'))
        try
            d = csvread(FileList{e});
            fidx = strfind(nm,'_');
            NM{end+1} = nm(1:fidx(end)-1);

            fileName = [p filesep NM{end} '_petioleData.csv'];
            K = csvread(fileName);
            
            fileName = [p filesep NM{end} '_rootProfile.csv'];
            R = csvread(fileName);
            
            D = [D;d];
            D1 = [D1;K];
            
            TOP = [TOP;d];
            BOT = [BOT;R'];
            
        catch
        end
    end
end
TOP = fliplr(TOP);
np = 1000;
for e = 1:size(BOT,1)
    [norLenBOT(e,:) nTotBOT(e,:)] = normalizeProfile(BOT(e,:),np);
    %nT(e,:) = normalizeProfile(TOP(e,:),np);
end
%% write out
for e = 1:numel(NM)
    forWIDTH{e,1} = NM{e};
    forWIDTH{e,2} = D1(e,3);
end
cell2csv('/mnt/spaldingdata/nate/communications/papers/carrotPaper/forWIDTH.csv', forWIDTH)
%% 2: load the 100 : step one: run before loading the 100
close all
[S C U E L ERR LAM] = PCA_FIT_FULL(TOP,5);
plot(C(:,1),C(:,2),'.')
cp1 = [1.5*10^4 0];
cp2 = [-1*10^4 2.2*10^4];
X = linspace(cp1(1),cp2(1),10);
Y = linspace(cp1(2),cp2(2),10);
hold on
plot(X,Y,'r*')
M = PCA_BKPROJ([X' Y'],E,U);
figure;plot(M')

%% 3: load the 100 : corr for 100
[w] = readtext('/home/nate/Downloads/petiole_hand_measurements.csv');
% get the header
HEADER = w(1,:);
w(1,:) = [];
NM2 = {};
leafNumber = [];
Y2 = [];
X = [];
N = [];
CT = [];
Co = [];
W = [];
% w is the hand measurement data
% [petCount newPetCount petWidth petW mean(petL) petNUM petDIA]
for e = 1:size(w,1)
    NM2{e} = w{e,1}(1:end-4);
    fidx =  find(strcmp(NM,NM2{e}));
    if ~isempty(fidx)
        % get leaf numerb
        leafNumber = [leafNumber;w{e,4}];
        % has pet length and width
        Y2 = [Y2;[w{e,2:3}]];
        % get lookup for 1 of 100 for TOP which is in D
        X = [X;D(fidx,:)];
        % get the pet width contained in D1
        N = [N;D1(fidx,3)];
        %
        CT = [CT;D1(fidx,[1 2 6])];
        Co = [Co;C(fidx,:)];
        % get all pet width
        W = [W;D1(fidx,[3 4 7])]
    end
    e
end
HEADER
%X = bsxfun(@times,X,N.^-1);
X100 = X;
%X = bsxfun(@times,X,Y2(:,2).^-1);
corr(leafNumber,CT);
corr(Y2,Co);
fprintf(['Measure all pet width profiles and select the best:\n'])
% Y2 : = hand measurements for pet width
% W := algorithm measurements for pet width
corr(Y2,W)
fprintf(['Measure all pet width profiles and select the best for hand:\n'])
corr(leafNumber,W)
%% 4: read manual pet counts 1000
% this will read in the whole data set for pet count/leaf count for most
% sidx is the filter for the data
% R is the profile data
% D1 comes from above and is the pet wdith - guess
petC = readtext('/home/nate/Downloads/manual_pet_counts.csv');
cnt = 1;
NM2 = {};
sidx = [];
for e = 2:size(petC,1)
    if ~ischar(petC{e,2})
        NM2{cnt} = petC{e,1}(1:(end-4));
        pv(cnt) = petC{e,2};
        cnt = cnt + 1;
    end
end
sidx = [];
R = [];
RV = [];
N = [];
for e = 1:numel(NM)
    fidx = find(strcmp(NM{e},NM2));
    if ~isempty(fidx)
        sidx(e) = fidx(1);
        RV = [RV;pv(fidx(1))];
        R = [R;TOP(e,:)];
        %R = [R;nT(e,:)];
        N = [N;D1(e,3)];
    else
        sidx(e) = 0;
    end
end
%% 5: load into REP data
XREP = bsxfun(@times,R,N.^-1);
YREP = RV;
%% 6: build pls regression from hand measurement PET COUNT
% note: for this to make sense - you must run the above and think
close all
COR = [];
Yp = [];
%{
% i think i did this wrong from before
% making wirht by changing the diviosn to the X
X = R;
Y = RV.*N.^-1; % divides the pet integration by the width chooose from the above

%{
STORE DATA FOR REPEAT
XREP = bsxfun(@times,R,N.^-1);
N = [ones(size(R,1)) R]*betaWIDTH];
XREP = bsxfun(@times,R,N.^-1);
YREP = RV;
%}
%}
X = XREP;
Y = YREP; % divides the pet integration by the width chooose from the above

[~ ,Co] = PCA_FIT_FULL(R,4);
%bLAM = diag(bLAM);
for l = 1:10
    parfor e = 1:size(X,1)
        idx = setdiff(1:size(X,1),e);
        
        subX = X(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        size(R);
        subY = Y(idx);
        
        %{
        %subX = [Co(idx,:) LV(idx,1)];
        subX = [Co(idx,:)];
        net = feedforwardnet([l 5]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        Yp(l,e) = sim(net,Co(e,:)');
        %}
        
        
        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
        %Yp = [ones(size(X,1),1) X(e,:)]*beta;
        %Yp(e) = [1 X(e,:)]*beta;
        %Yp(l,e) = [1 X(e,:) LV(e,1)]*beta;
        Yp(l,e) = [1 X(e,:)]*beta;
        
        %plot(Y(e),Yp(l,e),'.')
        drawnow
        hold on
    end
    COR(l) = corr(Y,Yp(l,:)','type','Pearson');
    l
    plot(COR)
    drawnow
end
close all
plot(COR)
[Xloadings,Yloadings,Xscores,Yscores,betaCOUNT,pctVar,mse,stats,Weights] = plsregress(X,Y,9);

%% 7: regress the pet count values
% from the ocean below
petCountRegressValue1000 = [ones(size(X,1),1) X]*betaCOUNT;
petCountRegressValue100 = [ones(size(X100,1),1) X100]*betaCOUNT;
%{
%% 8: build pls regression from hand measurement PET WIDTH - based on pet count
% note: for this to make sense - you must run the above and think
close all
COR = [];
Yp = [];
%{
% i think i did this wrong from before
% making wirht by changing the diviosn to the X
X = R;
Y = RV.*N.^-1; % divides the pet integration by the width chooose from the above
%}
X = bsxfun(@times,X100,petCountRegressValue100.^-1);
%X = X100;
Y = Y2(:,2);

[~ ,Co] = PCA_FIT_FULL(R,4);
%bLAM = diag(bLAM);
for l = 1:20
    parfor e = 1:size(X,1)
        idx = setdiff(1:size(X,1),e);
        
        subX = X(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        size(R);
        subY = Y(idx);
        
        %{
        %subX = [Co(idx,:) LV(idx,1)];
        subX = [Co(idx,:)];
        net = feedforwardnet([l 5]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        Yp(l,e) = sim(net,Co(e,:)');
        %}
        
        
        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
        %Yp = [ones(size(X,1),1) X(e,:)]*beta;
        %Yp(e) = [1 X(e,:)]*beta;
        %Yp(l,e) = [1 X(e,:) LV(e,1)]*beta;
        Yp(l,e) = [1 X(e,:)]*beta;
        
        %plot(Y(e),Yp(l,e),'.')
        drawnow
        hold on
    end
    COR(l) = corr(Y,Yp(l,:)');
    l
    plot(COR)
    drawnow
end
close all
[Xloadings,Yloadings,Xscores,Yscores,betaWIDTH,pctVar,mse,stats,Weights] = plsregress(X,Y,9);
%csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/pls_CORR_WID_over_model.csv',COR);
%WIDpre = [ones(size(X,1),1) X]*betaWIDTH;
%csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/pls_WIDSCATTER.csv',[Y WIDpre]);
%}
%% build pls regression for pet length
Y = Y2(:,1);
Yp = [];
close all
COR = [];
for l = 1:10
    parfor e = 1:size(X,1)
        idx = setdiff(1:size(X,1),e);
        %{
        subX = X(idx,:);
        %subX = Co(idx,:);
        %subX = nX(idx,:);
        subY = Y(idx);
        %}
        
        
        subY = Y(idx);
        subX = Co(idx,1:3);
        net = feedforwardnet([l 5]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        Yp(l,e) = sim(net,Co(e,1:3)');
        %Yp(e) = sim(net,X(e,:)');
        %
        
        
        %{
        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
        %Yp = [ones(size(X,1),1) X(e,:)]*beta;
        Yp(l,e) = [1 X(e,:)]*beta;
        %}
        %Yp(e) = [1 Co(e,:)]*beta;
        %plot(Y(e),Yp(e),'.')
        %drawnow
        %hold on
    end
    COR(l) = corr(Y,Yp(l,:)');
    plot(COR);
    drawnow
end
close all
plot(COR)
netLen = feedforwardnet([4 5]);
netLen.trainParam.showCommandLine = 1;
netLen.trainParam.showWindow = false;
netLen = train(netLen,Co(:,1:3)',Y');
lenPRE = sim(netLen,Co(:,1:3)');
plot(Y,lenPRE,'.')
csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/nn_petLength_predictions.csv',[Y lenPRE']);
csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/nn_petLength_predictions_CORR.csv',COR);
mean(abs(lenPRE-Y))
%% predict the 1000
lenPRE1000 = sim(netLen,C(:,1:3)');
for e = 1:numel(NM)
    forLENGTH{e,1} = NM{e};
    forLENGTH{e,2} = lenPRE1000(e);
end
cell2csv('/mnt/spaldingdata/nate/communications/papers/carrotPaper/forLENGTH.csv', forLENGTH)



%% PCA on data
[bS bC bU bE bL bERR bLAM] = PCA_FIT_FULL(BOT,4);
bLAM = diag(bLAM);
[nbS nbC nbU nbE nbL nbERR nbLAM] = PCA_FIT_FULL(norLenBOT,4);
nbLAM = diag(nbLAM);
[tnbS tnbC tnbU tnbE tnbL tnbERR tnbLAM] = PCA_FIT_FULL(nTotBOT,1000);
tnbLAM = diag(tnbLAM);
%% un normalized
close all
nC = 1;
[sweepD] = sweepPCA(bC,bE,bU,bLAM(1:nC).^.5,1:nC,5)
for e = 1:size(sweepD,1)
    figure
    for s = 1:size(sweepD,2)
        d = (squeeze(sweepD(e,s,:)));
        fidx = find(d < 2);
        d(fidx(1):end) = [];
        plot(d,-(1:numel(d)),'r');
        hold on
        plot(-d,-(1:numel(d)),'r');
    end
    axis equal
end
%% nor lenght
close all
nC = 2;
[sweepD] = sweepPCA(nbC,nbE,nbU,nbLAM(1:nC).^.5,1:nC,5);
for e = 1:size(sweepD,1)
    figure
    for s = 1:size(sweepD,2)
        d = (squeeze(sweepD(e,s,:)));
        %{
        fidx = find(d < 2);
        if ~isempty(fidx)
            d(fidx(1):end) = [];
        end
        %}
        plot(d,-(1:numel(d)),'r');
        hold on
        plot(-d,-(1:numel(d)),'r');
    end
    axis equal
end
%% nor lenght and width
close all
nC = 3;
mag = mean(max(BOT,[],1));
[sweepD] = sweepPCA(tnbC,tnbE,tnbU,1.5*tnbLAM(1:nC).^.5,1:nC,5);
for e = 1:size(sweepD,1)
    figure
    SWD = [];
    for s = 1:size(sweepD,2)
        d = (squeeze(sweepD(e,s,:)));
        %{
        fidx = find(d < 2);
        if ~isempty(fidx)
            d(fidx(1):end) = [];
        end
        %}
        d = mag*d;
        plot(d,-(1:numel(d)),'r');
        hold on
        plot(-d,-(1:numel(d)),'r');
        
        SWD = [SWD d];
    end
    csvwrite(['/mnt/spaldingdata/nate/communications/papers/carrotPaper/ROTT_PCASWEEP_' num2str(e) '.csv'], SWD)
    axis equal
end
%% top pre bot
close all
COR = [];
Yp = [];
[S C U E L ERR LAM] = PCA_FIT_FULL(BOT,4);
X = TOP;
Y = C;
comp = 2;
for l = 1:10
    parfor e = 1:size(X,1)
        idx = setdiff(1:size(X,1),e);
        %{
        subX = X(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        subY = Y(idx,comp);
        %}
        
        subX = [Co(idx,:)];
        net = feedforwardnet([l 5]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        Yp(l,e) = sim(net,[Co(e,:)]');
        
        %{
        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
        %Yp = [ones(size(X,1),1) X(e,:)]*beta;
        Yp(l,e) = [1 X(e,:)]*beta;
        %Yp(l,e) = [1 X(e,:) LV(e,1)]*beta;
        %}
        %plot(Y(e),Yp(e),'.')
        %drawnow
        hold on
    end
    COR(l) = corr(Y(:,comp),Yp(l,:)','type','Spearman');
end
close all
plot(COR)
%%
PL = [ones(size(TOP,1),1) TOP]*beta;
v
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(TOP,4);
[nbS nbC nbU nbE nbL nbERR nbLAM] = PCA_FIT_FULL(nB,4);
net = feedforwardnet([3 10 3]);
%net.trainParam.showCommandLine = 1;
%net.trainParam.showWindow = false;
%net = train(net,[C';PL'],bC');
net = train(net,bC',tC');
Yp = sim(net,bC');
%corr(Yp',bC)
M = PCA_BKPROJ(Yp',tE,tU);
close all

plot(M(10,:),'r')
hold on
plot(TOP(10,:),'k')
%% normalize profile
np = 1000;
for e = 1:size(BOT,1)
    nB(e,:) = normalizeProfile(BOT(e,:),np);
    nT(e,:) = normalizeProfile(TOP(e,:),np);
end
%% normalize the top
for e = 1:size(TOP,1)
    nT(e,:) = normalizeProfile(TOP(e,63:end),np);
end

%% from bot to top
close all
COR = [];
Yp = [];
[S C U E L ERR LAM] = PCA_FIT_FULL(TOP,4);
X = BOT;
Y = C;
comp = 2;
for l = 1:10
    parfor e = 1:size(X,1)
        idx = setdiff(1:size(X,1),e);
        
        subX = X(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        subY = Y(idx,comp);
        
        %{
        subX = [Co(idx,:) LV(idx,1)];
        net = feedforwardnet([l 5]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        Yp(l,e) = sim(net,[Co(e,:) LV(e,1)]');
        %}
        
        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(subX,subY,l);
        %Yp = [ones(size(X,1),1) X(e,:)]*beta;
        Yp(l,e) = [1 X(e,:)]*beta;
        %Yp(l,e) = [1 X(e,:) LV(e,1)]*beta;
        
        %plot(Y(e),Yp(e),'.')
        %drawnow
        hold on
    end
    COR(l) = corr(Y(:,comp),Yp(l,:)','type','Spearman');
end
close all
plot(COR)
%% get height from profile
for e = 1:size(TOP,1)
    he(e) = getProfileHeight(TOP(e,:))';
end
for e = 1:numel(he)
    HNM{e,1} = NM{e};
    HNM{e,2} = he(e);
end
cell2csv('/mnt/spaldingdata/nate/communications/papers/carrotPaper/HEIGHT.csv', HNM)
%% get wdith length from profile
heBOT = [];
for e = 1:size(BOT,1)
    heBOT(e,1) = getProfileHeight(BOT(e,:));
    heBOT(e,2) = max(BOT(e,:));
end
for e = 1:size(heBOT,1)
    HNMB{e,1} = NM{e};
    HNMB{e,2} = heBOT(e,1);
    HNMB{e,3} = heBOT(e,2);
end
cell2csv('/mnt/spaldingdata/nate/communications/papers/carrotPaper/ROOT_BLUNT.csv', HNMB)

%% double edge sword
X1 = R;
Y1 = RV;
X2 = X;
Y2 = W(:,1);
WIDTH = N;
WIDTH2 = W(:,1);
close all
iPLS2(X1,Y1,X2,Y2,4,4,WIDTH,WIDTH2);

%%
csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/pls_petCount_predictions.csv',Yp);
%% 
BB  = im2col(X,[1 50],'sliding');
[Sb Cb Ub Eb Lb ERRb LAMb] = PCA_FIT_FULL_T(BB,3);
for e = 1:size(Cb,1)
    nX(:,:,e) = col2im(Cb(e,:),[1 50],size(X));
end
sz = size(nX);
nX = reshape(nX,[sz(1) prod(sz(2:3))]);

%%
csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/nn_petLength_predictions.csv',Yp);
csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/pls_petLength_predictions.csv',Yp);
LV = Yp(4,:)';
%%
csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/PLS_CORR_pet_number.csv',COR);
csvwrite('/mnt/spaldingdata/nate/communications/papers/carrotPaper/PLS_COUNT_PRE.csv',[RV yPRE]);



