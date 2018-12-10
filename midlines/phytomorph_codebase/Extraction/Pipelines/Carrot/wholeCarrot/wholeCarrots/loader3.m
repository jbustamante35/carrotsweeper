%%%%%%%%%%%%%%%%%%%%%%
% 1: run this to load data from cyverse
%%%%%%%%%%%%%%%%%%%%%
FilePath = '/home/nate/Downloads/carrotModelBuild/wholeData/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
NM = {};
TOP = [];
BOT = [];
HIS = [];
CLASSIC_W = [];
CLASSIC_B = [];
CLASSIC_T = [];
% [petCount newPetCount petWidth petW mean(petL) petNUM petDIA]
for e = 1:numel(FileList)
    tic
    [p,nm,ext] = fileparts(FileList{e});
    
    
    if ~isempty(strfind(nm,'topProfile'))
        try
            %tmpTopProfile = csvread(FileList{e});
            fidx = strfind(nm,'_');
            NM{end+1} = nm(1:fidx(end)-1);
            
            fileName = [p filesep NM{end} '_whole_Perimeter.csv'];
            tmpWholeClassic(1) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_whole_Perimeter.csv'];
            tmpWholeClassic(1) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_whole_Solidity.csv'];
            tmpWholeClassic(2) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_whole_EulerNumber.csv'];
            tmpWholeClassic(3) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_whole_EquivDiameter.csv'];
            tmpWholeClassic(4) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_whole_Eccentricity.csv'];
            tmpWholeClassic(5) = csvread(fileName);
            
            CLASSIC_W = [CLASSIC_W;tmpWholeClassic];
            
            fileName = [p filesep NM{end} '_top_Perimeter.csv'];
            tmpy = csvread(fileName);
            tmptopClassic(1) = max(tmpy(:));
            
            fileName = [p filesep NM{end} '_top_Solidity.csv'];
            tmpy = csvread(fileName);
            tmptopClassic(2) = max(tmpy(:));
            
            fileName = [p filesep NM{end} '_top_EulerNumber.csv'];
            tmpy = csvread(fileName);
            tmptopClassic(3) = max(tmpy(:));
            
            fileName = [p filesep NM{end} '_top_EquivDiameter.csv'];
            tmpy = csvread(fileName);
            tmptopClassic(4) = max(tmpy(:));
            
            fileName = [p filesep NM{end} '_top_Eccentricity.csv'];
            tmpy = csvread(fileName);
            tmptopClassic(5) = max(tmpy(:));
            
            CLASSIC_T = [CLASSIC_T;tmptopClassic];
            
            fileName = [p filesep NM{end} '_bottom_Perimeter.csv'];
            tmpbottomClassic(1) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_bottom_Solidity.csv'];
            tmpbottomClassic(2) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_bottom_EulerNumber.csv'];
            tmpbottomClassic(3) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_bottom_EquivDiameter.csv'];
            tmpbottomClassic(4) = csvread(fileName);
            
            fileName = [p filesep NM{end} '_bottom_Eccentricity.csv'];
            tmpbottomClassic(5) = csvread(fileName);

            CLASSIC_B = [CLASSIC_B;tmpbottomClassic];
            
            
            % read the top profile
            fileName = [p filesep NM{end} '_topProfile.csv'];
            tmpTopProfile = csvread(fileName);
            % read the bottom profile
            fileName = [p filesep NM{end} '_rootProfile.csv'];
            tmpBottomProfile = csvread(fileName);
            
            
            fileName = [p filesep NM{end} '_histogramDistranceTransform.csv'];
            tmpH = csvread(fileName);
            
            HIS = [HIS;tmpH];
            
            TOP = [TOP;tmpTopProfile];
            BOT = [BOT;tmpBottomProfile];
            
        catch ME
            ME
            FileList{e}
            e
            break
        end
        toc
    end
end
%%%%%%%%%%%%%%%%%%%%%
%% 2: load the 100 : corr for 100
%%%%%%%%%%%%%%%%%%%%%
[w] = readtext('/home/nate/Downloads/petiole_hand_measurements.csv');
% get the header
HEADER = w(1,:);
w(1,:) = [];
NM2 = {};
leafNumber = [];
PET_WL = [];
histogramX100 = [];
profileTOPX100 = [];
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
        % get leaf nummber
        leafNumber = [leafNumber;w{e,4}];
        % has pet length and width
        PET_WL = [PET_WL;[w{e,2:3}]];
        % get lookup for 1 of 100 for TOP which is in D
        histogramX100 = [histogramX100;HIS(fidx,:)];
        profileTOPX100 = [profileTOPX100;TOP(fidx,:)];
        %{
        % get the pet width contained in D1
        N = [N;D1(fidx,3)];
        %
        CT = [CT;D1(fidx,[1 2 6])];
        Co = [Co;C(fidx,:)];
        % get all pet width
        W = [W;D1(fidx,[3 4 7])];
        %}
    end
    e
end
%%%%%%%%%%%%%%%%%%%%%
%% 3: read manual pet counts 1000
%%%%%%%%%%%%%%%%%%%%%
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

histogramX1000 = [];
profileTOPX1000 = [];
leafCount1000 = [];
%N = [];
for e = 1:numel(NM)
    fidx = find(strcmp(NM{e},NM2));
    if ~isempty(fidx)
        %if fidx(1) <= size(HIS,1)
            sidx(e) = fidx(1);

            leafCount1000 = [leafCount1000;pv(fidx(1))];
            %R = [R;TOP(e,:)];
            % get lookup for 1 of 100 for TOP which is in D
            histogramX1000 = [histogramX1000;HIS(e,:)];
            profileTOPX1000 = [profileTOPX1000;TOP(e,:)];
            %R = [R;nT(e,:)];
            %N = [N;D1(e,3)];
        %end
    else
        sidx(e) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%
%% 4: regress for width
%%%%%%%%%%%%%%%%%%%%%
h1 = figure;
h2 = figure;
rX = histogramX100;
%rX = [profileTOPX100 histogramX];
rY = PET_WL(:,2);
%[S rX U E L ERR LAM] = PCA_FIT_FULL(rX,5);
COR = [];
Yp = [];
for l = 1:10
    parfor e = 1:size(histogramX100,1)
        idx = setdiff(1:size(rX,1),e);
        
        subX = rX(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        %size(R);
        subY = rY(idx);
        
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
        Yp(l,e) = [1 rX(e,:)]*beta;
        %{
        figure(h2)
        plot(rY(e),Yp(l,e),'.')
        drawnow
        
        hold on
        %}
    end
    COR(l) = corr(rY(:),Yp(l,:)','type','Pearson');
    %l
    figure(h1)
    plot(COR)
    drawnow
end
[Xloadings,Yloadings,Xscores,Yscores,betaWIDTH,pctVar,mse,stats,Weights] = plsregress(rX,rY,2);
widthPrediction100 = [ones(size(histogramX100,1),1) histogramX100]*betaWIDTH;
widthPrediction1000 = [ones(size(histogramX1000,1),1) histogramX1000]*betaWIDTH;
close all
plot(rY(:),widthPrediction100,'.')
%%%%%%%%%%%%%%%%%%%%%
%% 5: regress for width - my regress in comment section
%%%%%%%%%%%%%%%%%%%%%
rX = histogramX100;
%rX = [profileTOPX100 histogramX];
rY = PET_WL(:,2);
%{
[S rX U E L ERR LAM] = PCA_FIT_FULL(rX,3);
P = eye(size(rX,2));
P = P(:);
func = @(P)findExpRegressorFunction(P,rX,rY);
options = optimoptions('fmincon','UseParallel',true,'Display','iter');
nP1 = fmincon(func,P,[],[],[],[],-2*ones(1,numel(P)),2*ones(1,numel(P)),[],options);
[cValue,betaWidth] = findExpRegressorFunction(nP1,rX,rY);
[rX] = PCA_REPROJ(histogramX100,E,U);
[Yp] = myNonPre(rX,nP1,betaWidth);
close all
plot(rY,real(Yp),'.')
corr(rY,real(Yp))
[C1000] = PCA_REPROJ(histogramX1000,E,U);
[widthPrediction10002] = myNonPre(C1000,nP1,betaWidth);
widthPrediction10002 = real(widthPrediction10002);


[C100] = PCA_REPROJ(histogramX100,E,U);
[widthPrediction1002] = myNonPre(C100,nP1,betaWidth);
widthPrediction1002 = real(widthPrediction1002);
%% regress for count
rX = profileTOPX1000;
%rX = [profileTOPX100 histogramX];
rY = leafCount1000(:);
[S rX U E L ERR LAM] = PCA_FIT_FULL(rX,3);
rX = [rX widthPrediction10002];

ridx = (rY > 20) | (widthPrediction10002 < 0);
rX(ridx,:) = [];
rY(ridx,:) = [];
rX = rX(1:200,:);
rY = rY(1:200);


P = eye(size(rX,2));
P = P(:);
func = @(P)findExpRegressorFunction(P,rX,rY);
options = optimoptions('fmincon','UseParallel',true,'Display','iter');
nP2 = fmincon(func,P,[],[],[],[],-2*ones(1,numel(P)),2*ones(1,numel(P)),[],options);
[cValue,betaCount] = findExpRegressorFunction(nP2,rX,rY);


[C2000] = PCA_REPROJ(profileTOPX1000,E,U);
rX = [C2000 widthPrediction10002];

rY = leafCount1000(:);
ridx = (rY > 20) | (widthPrediction10002 < 0);


rX(ridx,:) = [];
rY(ridx,:) = [];
rX = rX(201:end,:);
rY = rY(201:end);

[countPrediction10002] = myNonPre(rX,nP2,betaCount);
countPrediction10002 = real(countPrediction10002);
corr(countPrediction10002,rY)
close all

plot(countPrediction10002,rY,'.')
%}
%%%%%%%%%%%%%%%%%%%%%
%% 6: regress for count
%%%%%%%%%%%%%%%%%%%%%
close all
h1 = figure;
h2 = figure;
rX = bsxfun(@times,profileTOPX1000,widthPrediction1000.^-1);
%rX = profileTOPX1000;
rY = leafCount1000(:);
ridx = (rY > 30) | (widthPrediction1000 < 0);
rY(ridx) = [];
rX(ridx,:) = [];
COR = [];
Yp = [];
for l = 1:10
    parfor e = 1:size(rX,1)
        idx = setdiff(1:size(rX,1),e);
        
        subX = rX(idx,:);
        %subX = [X(idx,:) LV(idx,1)];
        %size(R);
        subY = rY(idx);
        
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
        Yp(l,e) = [1 rX(e,:)]*beta;
        %{
        figure(h2)
        plot(rY(e),Yp(l,e),'.')
        drawnow
        hold on
        %}
    end
    COR(l) = corr(rY(:),Yp(l,:)','type','Pearson');
    %l
    figure(h1)
    plot(COR)
    drawnow
end
[Xloadings,Yloadings,Xscores,Yscores,betaCOUNT,pctVar,mse,stats,Weights] = plsregress(rX,rY,9);
countPrediction1000 = [ones(size(rX,1),1) rX]*betaCOUNT;
close all
plot(rY(:),countPrediction1000,'.')
%%%%%%%%%%%%%%%%%%%%%
%% regress for length
%%%%%%%%%%%%%%%%%%%%%
h1 = figure;
h2 = figure;
[S C nn.U n.E L ERR LAM] = PCA_FIT_FULL(profileTOPX100,5);
rX = C;
rY = PET_WL(:,1);
COR = [];
Yp = [];
%rX = [profileTOPX100 histogramX];
%[S rX U E L ERR LAM] = PCA_FIT_FULL(rX,5);
for l = 1:10
    parfor e = 1:size(rX,1)
        idx = setdiff(1:size(rX,1),e);
        
        
        subY = rY(idx);
        subX = rX(idx,1:3);
        net = feedforwardnet([3 8]);
        net.trainParam.showCommandLine = 1;
        net.trainParam.showWindow = false;
        net = train(net,subX',subY');
        Yp(l,e) = sim(net,rX(e,1:3)');
        
        %{
        figure(h2)
        plot(rY(e),Yp(l,e),'.')
        drawnow
        hold on
        %}
    end
    COR(l) = corr(rY(:),Yp(l,:)','type','Pearson');
    %l
    figure(h1)
    plot(COR)
    drawnow
end
net = feedforwardnet([3 8]);
net.trainParam.showCommandLine = 1;
net.trainParam.showWindow = false;
net = train(net,rX',rY');
%%
genFunction(net,'/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Carrot/wholeCarrot/wholeCarrots/petLength.m');
nn.func = @(X)petLength(X);
%%
stage2Func = @(X)singleWholeCarrotStage2(X,betaWIDTH,betaCOUNT,nn);
pF = partialFunction(stage2Func,'petLength');
pF.publish();
%%
 websave('/home/nate/Downloads/petLength.mat','https://de.cyverse.org/dl/d/615CA1A4-CD2F-4F25-9E8E-0AA0BDCD9639/petLength.mat');
                    load('/home/nate/Downloads/petLength.mat');
                    func = obj.func;


