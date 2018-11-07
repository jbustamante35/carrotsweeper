%% write out csv data again
FilePath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/new_return/mat/';
FilePath = '/mnt/piranhapuffer2/Hannah/Maize QTL Gravi/return/mat/';
FileList = {};
FileExt = {'mat'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
% init vars for data
SNIP = 20;
NP = 20;
NPK = 20;
toWrite = 1;
outPath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/new_return/';
outPath = '/mnt/piranhapuffer2/Hannah/Maize QTL Gravi/return/';

for e = 1:numel(FileList)
    % save to mat file
    [pth,nm,ext] = fileparts(FileList{e});
    pth = nm;
    if toWrite
        % write out data
        csvOutPath = [outPath 'csv2/'];
        mkdir(csvOutPath);
        writeData(FileList{e},[csvOutPath pth],NP,SNIP,NPK);
    end
end
%% dig for csv
FilePath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/new_return/csv/';
FilePath = '/mnt/spaldingdata/nate/20-Nov-2013_14_58_27/';
FilePath = '/mnt/piranhapuffer2/Hannah/Maize QTL Gravi/return/csv/';
FileList = {};
FileExt = {'csv'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% co-plot data from WT compare
FR = 240;
clear data
close all
sam = 30;
searchStr = {};
searchStr{1} = {'MQW1','MQW2','MQW3','MQW4','MQW5','MQW6','MQW7'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 1;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'Wild Type';
data(n).color = 'k';
data(n).lineWidth = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ1'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 2;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'MQ1';
data(n).color = 'b';
data(n).lineWidth = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ2'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 3;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'MQ2';
data(n).color = 'r';
data(n).lineWidth = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ3'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 4;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'MQ3';
data(n).color = 'g';
data(n).lineWidth = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ4'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 5;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'MQ4';
data(n).color = 'c';
data(n).lineWidth = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ5'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 6;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'MQ5';
data(n).color = 'b';
data(n).lineWidth = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ6'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 7;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'MQ6';
data(n).color = 'r';
data(n).lineWidth = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ7'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = 8;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'MQ7';
data(n).color = 'g';
data(n).lineWidth = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'};
searchStr{2} = {'angle'};
searchStr{3} = {'10June14','11June14','16June14','17June14','18June14','22June14','23June14'};
n = numel(data)+1;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'all';
data(n).color = 'r';
data(n).lineWidth = 3;

mainF = figure;
for e =1:numel(data)
    data(e).U = mean(data(e).data,2);
    data(e).STD = std(data(e).data,1,2);
    data(e).SE = data(e).STD.*size(data(e).data,2).^-.5;
end
LEG = {};
for e = 1:numel(data)
    TIME = 1:size(data(e).data,1);
    errorbar(TIME(1:sam:end),data(e).U(1:sam:end),data(e).SE(1:sam:end),data(e).color,'LineWidth',data(e).lineWidth);
    hold all
    LEG{e} = [data(e).name '--' num2str(size(data(e).data,2))];
end

totD = [data.data];
labs = [data.label];
[S C U E L ERR LAM] = PCA_FIT_FULL(totD',8);
WTgn = 1;
UQ = unique(labs);
wtidx = find(UQ==WTgn);
UQ(wtidx) = [];
clear Rhv Rpv lambda lambda_vec
for u = 1:numel(UQ)
    mtg = labs == UQ(u);
    wtg = labs == WTgn;    
    tmpL = labs((wtg | mtg));
    tmpD = C((wtg | mtg),:);
    tmprD = totD(:,(wtg | mtg))';
    [tS tC tU{u} tE tL tERR tLAM] = PCA_FIT_FULL(tmprD,4);
    tmpD = tC;
    tmp = mynLDA(tmpD,tmpL,1,2);    
    lambda(:,u) = tmp(:,1);
    
    lambda_vec(u,:) = PCA_BKPROJ(lambda(:,u)',tE,0);
    num{u} = tmpD*lambda(:,u);
    
    [hv(u) pv(u)] = ttest2(num{u}(tmpL==WTgn,1),num{u}(tmpL==UQ(u),1));
    [Rhv(u,:) Rpv(u,:)] = ttest2(tmprD(tmpL==WTgn,:),tmprD(tmpL==UQ(u),:),[],[],[],1);
    plot(TIME,-Rhv(u,:)*u/10,[data(UQ(u)).color '-'],'LineWidth',3);    
    LEG{end+1} = ['raw ttest'];
end

for e = 1:30
    ulambda_vec = mean(lambda_vec,1);
    for u = 1:numel(UQ)
        if u ~= 1
            if sign(ulambda_vec*lambda_vec(u,:)') < 1
                lambda_vec(u,:) = -lambda_vec(u,:);
            end
        end       
    end
end

subF = figure;hold on;
nonBindVec = figure;hold on;
title('Orginal LDA vector');
ulambda_vec = mean(lambda_vec,1);
for u = 1:numel(UQ)
    % create labels
    mtg = labs == UQ(u);
    wtg = labs == WTgn;    
    tmpL = labs((wtg | mtg));
    NUM_WT_u = mean(num{u}(tmpL==WTgn));
    NUM_MT_u = mean(num{u}(tmpL==UQ(u)));
    
    WTCURVE = NUM_WT_u*lambda_vec(u,:) + tU{u};
    MTCURVE = NUM_MT_u*lambda_vec(u,:) + tU{u};
    
    
    
    figure(subF);
    %plot(WTCURVE,'k');
    %plot(MTCURVE,data(u+1).color);    
    plot(WTCURVE - MTCURVE,data(UQ(u)).color);
    
    figure(nonBindVec)
    plot(lambda_vec(u,:),data(UQ(u)).color);
    
    lambda_vecB(u,:) = bindVec(lambda_vec(u,:));    
    figure(mainF);
    plot(TIME,lambda_vecB(u,:)*pi/2,[data(UQ(u)).color '-'],'LineWidth',data(UQ(u)).lineWidth);
    LEG{end+1} = [num2str(hv(u)) '@' num2str(pv(u))];
end
legend(LEG);
%% generate bar chart from WT compare
figure;
mU = [];
mS = [];

for u = 1:numel(UQ)
     % create labels
    mtg = labs == UQ(u);
    wtg = labs == WTgn;    
    tmpL = labs((wtg | mtg));
    U = [mean(num{u}(tmpL==WTgn)) mean(num{u}(tmpL==UQ(u)))];
    S = [std(num{u}(tmpL==WTgn))*sum(tmpL==WTgn)^-.5 std(num{u}(tmpL==UQ(u)))*sum(tmpL==UQ(u))^-.5];    
    U = U - mean(U);
    mU = [mU U];
    mS = [mS S];
end


[hBar hErrorbar] = barwitherr(mS,1:numel(mU),mU,1);
for u = 1:numel(UQ)
    set(hBar(2*u),'FaceColor',data(UQ(u)).color);
    set(hBar(2*u-1),'FaceColor','w');
end
%% one to split them all
UQ = 1:numel(data);
clear proj
clear U S
for u = 1:numel(UQ)
    proj{u} = data(u).data'*lambda_vec(end,:)';
    U(u) = mean(proj{u});
    S(u) = std(proj{u})*numel(proj{u})^-.5;
end
U = U - mean(U);
figure;
barwitherr(S,U)
%% co-plot data from WS
FR = 240;
clear data
close all
sam = 1;
searchStr = {};
searchStr{1} = {'WS'};
searchStr{2} = {'angle'};
n = 1;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = 'Wild Type - WS';
data(n).color = 'k';
data(n).lineWidth = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
searchStr = {};
searchStr{1} = {'_2_1'};
searchStr{2} = {'angle'};
n = 2;
[data(n).data] = simpleLoader(searchStr,FileList,FR);
data(n).label = n*ones(1,size(data(n).data,2));
[data(n).data data(n).label] = simpleCleaner(data(n).data,data(n).label,15,40);
data(n).name = '2-1';
data(n).color = 'b';
data(n).lineWidth = 1;




mainF = figure;
for e =1:numel(data)
    data(e).U = mean(data(e).data,2);
    data(e).STD = std(data(e).data,1,2);
    data(e).SE = data(e).STD.*size(data(e).data,2).^-.5;
end
LEG = {};
for e = 1:numel(data)
    TIME = 1:size(data(e).data,1);
    errorbar(TIME(1:sam:end),data(e).U(1:sam:end),data(e).SE(1:sam:end),data(e).color,'LineWidth',data(e).lineWidth);
    hold all
    LEG{e} = [data(e).name '--' num2str(size(data(e).data,2))];
end

totD = [data.data];
labs = [data.label];
[S C U E L ERR LAM] = PCA_FIT_FULL(totD',8);
WTgn = 1;
UQ = unique(labs);
wtidx = find(UQ==WTgn);
UQ(wtidx) = [];
clear Rhv Rpv lambda lambda_vec
for u = 1:numel(UQ)
    mtg = labs == UQ(u);
    wtg = labs == WTgn;    
    tmpL = labs((wtg | mtg));
    tmpD = C((wtg | mtg),:);
    tmprD = totD(:,(wtg | mtg))';
    [tS tC tU{u} tE tL tERR tLAM] = PCA_FIT_FULL(tmprD,4);
    tmpD = tC;
    tmp = mynLDA(tmpD,tmpL,1,2);    
    lambda(:,u) = tmp(:,1);
    
    lambda_vec(u,:) = PCA_BKPROJ(lambda(:,u)',tE,0);
    num{u} = tmpD*lambda(:,u);
    
    [hv(u) pv(u)] = ttest2(num{u}(tmpL==WTgn,1),num{u}(tmpL==UQ(u),1));
    [Rhv(u,:) Rpv(u,:)] = ttest2(tmprD(tmpL==WTgn,:),tmprD(tmpL==UQ(u),:),[],[],[],1);
    plot(TIME,-Rhv(u,:)*u/10,[data(UQ(u)).color '-'],'LineWidth',3);    
    LEG{end+1} = ['raw ttest'];
end

for e = 1:30
    ulambda_vec = mean(lambda_vec,1);
    for u = 1:numel(UQ)
        if u ~= 1
            if sign(ulambda_vec*lambda_vec(u,:)') < 1
                lambda_vec(u,:) = -lambda_vec(u,:);
            end
        end       
    end
end

subF = figure;hold on;
nonBindVec = figure;hold on;
title('Orginal LDA vector');
ulambda_vec = mean(lambda_vec,1);
for u = 1:numel(UQ)
    % create labels
    mtg = labs == UQ(u);
    wtg = labs == WTgn;    
    tmpL = labs((wtg | mtg));
    NUM_WT_u = mean(num{u}(tmpL==WTgn));
    NUM_MT_u = mean(num{u}(tmpL==UQ(u)));
    
    WTCURVE = NUM_WT_u*lambda_vec(u,:) + tU{u};
    MTCURVE = NUM_MT_u*lambda_vec(u,:) + tU{u};
    
    
    
    figure(subF);
    %plot(WTCURVE,'k');
    %plot(MTCURVE,data(u+1).color);    
    plot(WTCURVE - MTCURVE,data(UQ(u)).color);
    
    figure(nonBindVec)
    plot(lambda_vec(u,:),data(UQ(u)).color);
    
    lambda_vecB(u,:) = bindVec(lambda_vec(u,:));    
    figure(mainF);
    plot(TIME,lambda_vecB(u,:)*pi/2,[data(UQ(u)).color '-'],'LineWidth',data(UQ(u)).lineWidth);
    LEG{end+1} = [num2str(hv(u)) '@' num2str(pv(u))];
end
legend(LEG);

%%
figure;
searchStr = {};
%searchStr{1} = 'vps2';
searchStr{1} = {'lip5-1'};
searchStr{1} = {'MQ1'};
searchStr{2} = {'length'};
%searchStr{3} = {'10June14'};
%searchStr{3} = '12';
%searchStr{4} = 'Dec';
ret = find(strcontains(FileList,searchStr));
MT = [];
FR = 241;
for e = 1:numel(ret)
    fn = FileList{ret(e)};
    D = csvread(fn);
    MT = [MT D(1:FR,:)];
end
%ridx = find(any(MT == 0,1) | any(abs(diff(MT,1,1)) > 10,1) | abs(MT(1,:)) > 40);
%MT(:,ridx) = [];
MT = bsxfun(@minus,MT,MT(1,:));
plot(MT)
%%
searchStr = {};
%searchStr{1} = 'vps2';
searchStr{1} = {'lip5-2'};
searchStr{2} = {'angle'};
searchStr{2} = {'length'};
%searchStr{3} = '12';
%searchStr{4} = 'Dec';
ret = find(strcontains(FileList,searchStr));
MT2 = [];
FR = 360;
for e = 1:numel(ret)
    fn = FileList{ret(e)};
    D = csvread(fn);
    MT2 = [MT2 D(1:FR)];
end
ridx = find(any(MT2 == 0,1) | any(abs(diff(MT2,1,1)) > 10,1) | abs(MT2(1,:)) > 40);
MT2(:,ridx) = [];
plot(MT2)
%%
close all
uWT = mean(WT,2);
sWT = std(WT,1,2)*size(WT,2)^-.5;
uMT = mean(MT,2);
sMT = std(MT,1,2)*size(MT,2)^-.5;
uMT2 = mean(MT2,2);
sMT2 = std(MT2,1,2)*size(MT2,2)^-.5;

[h p] = ttest2(WT,MT,[],[],[],2);
[h2 p] = ttest2(WT,MT2,[],[],[],2);
figure;
TIME = linspace(1,360*2/60,360);
errorbar(TIME,uWT,sWT,'k');
hold on
errorbar(TIME,uMT,sMT,'r');
errorbar(TIME,uMT2,sMT2,'b');
axis([1 TIME(end) 0 100]);
tt = axes;
plot(TIME,.4*h,'r');
hold on
plot(TIME,.8*h2,'b');
set(tt,'YAxisLocation','right','Color','none')
axis([1 TIME(end) 0 2]);
csvwrite('/mnt/spaldingdata/nate/20-Nov-2013_14_58_27/WT_LIP5-1.csv',[uWT sWT uMT sMT h]);
csvwrite('/mnt/spaldingdata/nate/20-Nov-2013_14_58_27/WT_LIP5-2.csv',[uWT sWT uMT2 sMT2 h2]);
%%

close all
uWT = mean(WT,2);
sWT = std(WT,1,2)*size(WT,2)^-.5;
uMT = mean(MT,2);
sMT = std(MT,1,2)*size(MT,2)^-.5;


[h p] = ttest2(WT,MT,[],[],[],2);

figure;
TIME = linspace(1,size(uWT,1)*2/60,size(uWT,1));
errorbar(TIME,uWT,sWT,'k');
hold on
errorbar(TIME,uMT,sMT,'r');
%axis([1 TIME(end) 0 100]);
%tt = axes;
plot(TIME,max(uWT(:))*h/2,'r');
hold on
set(tt,'YAxisLocation','right','Color','none')
axis([1 TIME(end) 0 2]);
%%
totD = [WT MT];
labs = [zeros(1,size(WT,2)) ones(1,size(MT,2))];
[S C U E L ERR LAM] = PCA_FIT_FULL(totD',4);
[lambda] = mynLDA(C,labs,1,3);
lambda_vec = PCA_BKPROJ(lambda',E,0);
lambda_vec = bindVec(lambda_vec);
plot(TIME,lambda_vec*max(uWT(:)));
num = C*lambda;
for g = 1:size(num,2)    
    [pv(g) hv(g)] = ttest2(num(find(labs==0),g),num(labs==1,g));
end
%% load K
FilePath = '/mnt/piranhapuffer2/Hannah/Maize QTL Gravi/return/mat/';
FileList = {};
FileExt = {'mat'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
% init vars for data
SNIP = 20;
NP = 20;
NPK = 150;

%%
clear out
parfor e = 1:numel(FileList)
    try
        tm = clock;
        out{e} = loadData(FileList{e},SNIP,NP,NPK);
        fprintf(['Done@' num2str(e) ':' num2str(etime(clock,tm)) '\n']);
    catch ME
        ME
    end        
end
%%
close all
searchStr = {};
searchStr{1} = {'MQW'};

searchStr2{1} = {'MQ3','MQ7'};
WTKUR = [];
MTKUR = [];
FR = 240;
for e = 1:numel(out)
    if strcontains({out{e}.dataID},searchStr)
        try
        WTKUR = cat(3,WTKUR,out{e}.K(:,1:FR,:));
        catch
        end
    end
    
    if strcontains({out{e}.dataID},searchStr2)
        try
            MTKUR = cat(3,MTKUR,out{e}.K(:,1:FR,:));
        catch
        end
    end
end
figure;
mesh(mean(WTKUR,3))
view([0 90])
figure;
mesh(mean(MTKUR,3))
view([0 90])
[P H] = ttest2(WTKUR,MTKUR,[],[],[],3);
figure;
mesh(H)
view([0 90]);
figure;
mesh(P)
view([0 90]);

