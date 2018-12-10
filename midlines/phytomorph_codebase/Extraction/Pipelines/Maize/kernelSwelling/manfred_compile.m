% read in data and create mastersheet
mPath = '/mnt/snapper/Manfred/return/';
%mPath = '/mnt/snapper/kernelSwellingData/Jeff/return/';
[MS rL] = createMasterFile(mPath);
outName = '/mnt/snapper/Manfred/masterSheetv2.csv';
% COL--1:2:para--3:error--4:291:swell---292:580:fit--581:869--area
cell2csv(outName,MS);
%% lookup key by folder
close all
clear fidx
strSWELL = 12:(12+288-1);
fidx{1} = find(strcmp(MS(:,1),'Mon') & strcmp(MS(:,2),'Jan') & strcmp(MS(:,3),'23') & strcmp(MS(:,4),'12-28-57'));
fidx{2} = find(strcmp(MS(:,1),'Mon') & strcmp(MS(:,2),'Jan') & strcmp(MS(:,3),'23') & strcmp(MS(:,4),'12-47-13'));
fidx{3} = find(strcmp(MS(:,1),'Mon') & strcmp(MS(:,2),'Jan') & strcmp(MS(:,3),'23') & strcmp(MS(:,4),'13-08-45'));
gidx = find(strcmp(MS(:,1),'Fri') & strcmp(MS(:,2),'Feb') & strcmp(MS(:,3),'10'));
sw1 = cell2mat(MS(fidx{1},strSWELL));
sw2 = cell2mat(MS(fidx{2},strSWELL));
sw3 = cell2mat(MS(fidx{3},strSWELL));
plot(mean(sw1,1),'r')
hold on
plot(mean(sw2,1),'g')
plot(mean(sw3,1),'b')
%% create look up tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:size(MS,1)
    keyT{e} = [MS{e,[1 2 3 5 6 7]}];
end

LOOKUPTABLE = fliplr(1:16);

GN = readtext('/home/nate/Downloads/Landraces_imbibition_warm_rep1.csv');
for e = 1:size(GN,1)
    for c = 1:size(GN,2)
        if ~isstr(GN{e,c})

            if c == 6
                GN{e,c} = LOOKUPTABLE(GN{e,c});
            end

            if GN{e,c} < 10 & c ~= 5
                pad = '0';
            else
                pad  = '';
            end
            GN{e,c} = [pad num2str(GN{e,c})];


        end
    end
    keyL{e} = [GN{e,1:end-1}];
    keySS{e} = [GN{e,1} '--' GN{e,5} '--' GN{e,6}];
    keyG{e} = GN{e,end};
end

for e = 1:numel(keyT)
    fidx = find(strcmp(keyT{e},keyL));
    keyGT{e} = keyG{fidx};
    keySST{e} = keySS{fidx};
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QC on contours
mPath = '/mnt/snapper/Manfred/return/';
mFileList = {};
FileExt = {'mat'};
verbose = 1;
mFileList = gdig(mPath,mFileList,FileExt,verbose);
%%
close all
sel = 280;
a = load(mFileList{sel});
fidx = strfind(mFileList{sel},'--');
BOX = 200;
SKIP = 5;
kernelNumber = 8;
dataPath = ['/mnt/snapper/Manfred/rawData/' mFileList{sel}((fidx(1)+2):fidx(2)-1) filesep mFileList{sel}((fidx(2)+2):fidx(3)-1) filesep];
for t = 1:SKIP:288
    imgFile = [dataPath num2str(t) '.tif'];
    tI = imread(imgFile,'PixelRegion',{[a.tmpCenters(kernelNumber,1)-BOX a.tmpCenters(kernelNumber,1)+BOX],[a.tmpCenters(kernelNumber,2)-BOX a.tmpCenters(kernelNumber,2)+BOX]});
    imshow(tI,[])
    tC = a.tmpbdB{kernelNumber}{t};
    if t == 1
        iC = tC;
    end
    %tC = bsxfun(@minus,tC,a.tmpCenters(kernelNumber,:)-BOX)
    hold on
    plot(tC(:,2),tC(:,1),'m');
    plot(iC(:,2),iC(:,1),'r')
    drawnow
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at all B73
strSWELL = 12;
strPARA = 9:10;
strAREA = 589:595;
tmp = cell2mat(MS(:,strSWELL:(strSWELL+288-1)));
para = cell2mat(MS(:,strPARA));
areaT = cell2mat(MS(:,strAREA));
for e = 1:size(MS,1)
    key{e} = [MS{e,1:7}];
end
toFit = 1;
clear U S
UQ = unique(key);
cnt = 1;
for u =1:numel(UQ)    
    fidx = strcmp(key,UQ{u});
    
    if strcmp(unique(keyGT(fidx)),'B73')

        
        
        sub = tmp(fidx,1:end);
        subPara = para(fidx,:);
        [keep] = manfredFilter(sub);
        
        if toFit
            sub = func(subPara(:,1),subPara(:,2),1:500);
        end
        
        plot(sub','k')
        hold on
        plot(sub(keep,:)','r')
        hold off
        axis([0 size(sub,2) 0 .2])
        %waitforbuttonpress
        U(cnt,:) = mean(sub(keep,:),1);
        S(cnt,:) = std(sub(keep,:),1,1)*sum(keep)^-.5;
        cnt = cnt + 1;
    end
end

close all
errorbar(U',S')



%% group by genotype and get swell rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = cell2mat(MS(:,strSWELL:(strSWELL+288-1)));
para = cell2mat(MS(:,strPARA));
areaT = mean(cell2mat(MS(:,strAREA)),2);
close all
clear U S P
UQ = unique(keyGT);
for u = 1:numel(UQ)
    fidx = strcmp(keyGT,UQ{u});
    
    sub = tmp(fidx,1:end);
    subPara = para(fidx,:);
    
    [keep] = manfredFilter(sub);

    if toFit
        sub = func(subPara(:,1),subPara(:,2),1:500);
    end
    subArea = areaT(fidx,:);
    A(u,:) = mean(subArea(keep));
    P(u,:) = mean(subPara(keep,:),1);    
    %{
    plot(sub','k')
    plot(sub(keep,:)','r')    
    axis([0 288 0 .2])
    %}
    
    U(u,:) = mean(sub(keep,:),1);
    S(u,:) = std(sub(keep,:),1,1)*sum(keep)^-.5;
    errorbar(U(u,:),S(u,:))
    hold all
    title(UQ{u})
    %drawnow
    %waitforbuttonpress
end
close all
errorbar(U',S')
%% sort by rate run above first
close all
fun = figure;
plot3(P(:,1),P(:,2),A,'k.')
hold on
old = figure;
[swellRate,sidx] = sort(P(:,2));
sortGenoByRate = UQ(sidx);
topN = 5;
slowest = sortGenoByRate(1:topN);
fastest = sortGenoByRate((end-(topN-1)):end);
mid = sortGenoByRate(148:152);
plot(U','k')
hold on
clear U S
toFit = 1;
for u = 1:numel(slowest)
    fidx = strcmp(keyGT,slowest{u});
    
    sub = tmp(fidx,1:end);
    subPara = para(fidx,:);
    [keep] = manfredFilter(sub);

    if toFit
        sub = func(subPara(:,1),subPara(:,2),1:500);
    end
    
    subArea = areaT(fidx,:);
    uA = mean(subArea(keep));
    
    uP = mean(subPara(keep,:),1);
    figure(fun)
    plot3(uP(1),uP(2),uA,'ro');
    %axis([0 1 0 .015])
    figure(old);    
    
    
    %{
    plot(sub','k')
    plot(sub(keep,:)','r')    
    axis([0 288 0 .2])
    %}
    
    U(u,:) = mean(sub(keep,:),1);
    S(u,:) = std(sub(keep,:),1,1)*sum(keep)^-.5;
    %errorbar(U(u,:),S(u,:),'r')
    plot(U(u,:),'r','LineWidth',3)
    hold all
    %plot(sub(keep,:)','r')
    hold all
    title(UQ{u})
    drawnow
    %waitforbuttonpress
end

for u = 1:numel(mid)
    fidx = strcmp(keyGT,mid{u});
    
    sub = tmp(fidx,1:end);
    subPara = para(fidx,:);
    [keep] = manfredFilter(sub);

    if toFit
        sub = func(subPara(:,1),subPara(:,2),1:500);
    end
    
    subArea = areaT(fidx,:);
    uA = mean(subArea(keep));
    
    uP = mean(subPara(keep,:),1);
    figure(fun)
    plot3(uP(1),uP(2),uA,'bo');
    figure(old);
    [keep] = manfredFilter(sub);
    %{
    plot(sub','k')
    plot(sub(keep,:)','r')    
    axis([0 288 0 .2])
    %}
    
    U(u,:) = mean(sub(keep,:),1);
    S(u,:) = std(sub(keep,:),1,1)*sum(keep)^-.5;
    %errorbar(U(u,:),S(u,:))
    plot(U(u,:),'b','LineWidth',3)
    %errorbar(U(u,:),S(u,:),'b')
    hold all
    title(UQ{u})
    drawnow
    %waitforbuttonpress
end

for u = 1:numel(fastest)
    fidx = strcmp(keyGT,fastest{u});
    
    sub = tmp(fidx,1:end);
    subPara = para(fidx,:);
    [keep] = manfredFilter(sub);

    if toFit
        sub = func(subPara(:,1),subPara(:,2),1:500);
    end
    
    subArea = areaT(fidx,:);
    uA = mean(subArea(keep));
    
    uP = mean(subPara(keep,:),1);
    figure(fun)
    %plot(uP(1),uP(2),'go');
    plot3(uP(1),uP(2),uA,'go');
    figure(old);
    [keep] = manfredFilter(sub);
    %{
    plot(sub','k')
    plot(sub(keep,:)','r')    
    axis([0 288 0 .2])
    %}
    
    U(u,:) = mean(sub(keep,:),1);
    S(u,:) = std(sub(keep,:),1,1)*sum(keep)^-.5;
    %errorbar(U(u,:),S(u,:))
    plot(U(u,:),'g','LineWidth',3)
    %errorbar(U(u,:),S(u,:),'g')
    hold all
    title(UQ{u})
    drawnow
    %waitforbuttonpress
end

plot(CONU,'b','LineWidth',3)
%% indiv
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fidx = strcmp(keyGT,'PE0131');
sub = tmp(fidx,1:end);
if toFit
    subPara = para(fidx,:);
    sub = func(subPara(:,1),subPara(:,2),1:500);
end
d_sub = abs(diff(sub,1,2));


keep = ~any(d_sub > .3,2);
CONU = mean(sub(keep,:),1);
CONS = std(sub(keep,:),1,1)*sum(keep)^-.5;



tmp = cell2mat(MS(:,strSWELL:(strSWELL+288-1)));
%close all
figure;
clear U S
UQ = unique(keyGT);
for u =1:numel(UQ)
    fidx = strcmp(keyGT,UQ{u});
    sub = tmp(fidx,1:end);
    subPara = para(fidx,:);
    [keep] = manfredFilter(sub);

    if toFit
        sub = func(subPara(:,1),subPara(:,2),1:500);
    end
    
    plot(sub','k')
    hold on
    plot(sub(keep,:)','r')    
    axis([0 size(sub,2) 0 .4])
    plot(CONU,'g','LineWidth',3)
    
    U(u,:) = mean(sub(keep,:),1);
    S(u,:) = std(sub(keep,:),1,1)*sum(keep)^-.5;
    %errorbar(U(u,:),S(u,:),'b')
    plot(U(u,:),'b')
    
    
    
    %title(UQ{u})
    REP = unique(keySST(fidx));
    title([UQ{u} '--' REP])
    %waitforbuttonpress
    hold off
    drawnow
end

%close all
figure;
errorbar(U',S')
figure;
plot(U','k')
hold on
plot(CONU,'g','LineWidth',3)
%% add by genoType list

close all
errorbar(U',S')
figure;
plot(U','k')
hold on
plot(CONU,'g','LineWidth',3)
gL = {'PE0029','PE0001','KE0173','PE0161','KE0286','KE0149'};
gL = {'PE0060','PE0070','KE0005','PE0029','PE0071','PE0023','B73'};
subS = [];
for g = 1:numel(gL)
    fidx = find(strcmp(gL{g},keyGT));
    sub = tmp(fidx,1:end);
    subPara = para(fidx,:);
    
    [keep] = manfredFilter(sub);

    
    if toFit
        sub = func(subPara(:,1),subPara(:,2),1:500);
    end
    
    subSR(g) = mean(subPara(keep,2));
    
    subS = [subS;subPara(keep,2)];
    
    
    plot(mean(sub(keep,:),1),'LineWidth',3)
    hold all
end


MS(:,1)
figure;
ksdensity(P(:,2));
hold on
ksdensity(subS);

tmpD = [];
kidx = find(strcmp(MS(:,1),'Fri'));
tmpK = setdiff(keyGT,gL);
for p = 1:numel(keyGT)
    if any(strcmp(keyGT{p},tmpK))
        if MS{p,10} < .04
            tmpD = [tmpD ;MS{p,10}];
        end
    end
end
figure;
%tmpD = cell2mat(MS(kidx,10));
%tmpD(tmpD > .04) = [];
ksdensity(tmpD)
hold on
ksdensity(subS);

uP = mean(P(:,2),1);
sP = std(P(:,2),1,1);
suP = mean(subSR',1);
ssP = std(subSR',1,1);
figure;
errorbar([uP suP],[sP ssP])

%% by land race
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fidx = strcmp(keyGT,'B73');
sub = tmp(fidx,5:end);
d_sub = abs(diff(sub,1,2));
keep = ~any(d_sub > .3,2);
CONU = mean(sub(keep,:),1);
CONS = std(sub(keep,:),1,1)*sum(keep)^-.5;


tmp = cell2mat(MS(:,strSWELL:(strSWELL+288-1)));
close all
clear U S
UQ = unique(keyGT);
for e = 1:3
    LR{e} = [];
end
for u =1:numel(UQ)
    fidx = strcmp(keyGT,UQ{u});
    
    sub = tmp(fidx,1:end);
    subPara = para(fidx,:);
    [keep] = manfredFilter(sub);

    if toFit
        sub = func(subPara(:,1),subPara(:,2),1:500);
    end
    
    plot(sub','k')
    hold on
    plot(sub(keep,:)','r')    
    axis([0 size(sub,2) 0 .4])
    plot(CONU,'g','LineWidth',3)
    
    U(u,:) = mean(sub(keep,:),1);
    S(u,:) = std(sub(keep,:),1,1)*sum(keep)^-.5;
    %errorbar(U(u,:),S(u,:),'b')
    plot(U(u,:),'b')
    
    switch UQ{u}(1:2)
        case 'KE'
            LR{1} = [LR{1};sub(keep,:)];
        case 'PE'
            LR{2} = [LR{2};sub(keep,:)];
        case 'LL'
            LR{3} = [LR{3};sub(keep,:)];
    end
    REP = unique(keySST(fidx));
    title([UQ{u} '--' REP])
    %waitforbuttonpress
    hold off
    drawnow
end
close all
errorbar(U',S')
figure;
plot(U','k')
hold on
plot(CONU,'g','LineWidth',3)
close all
CL = {'r' 'g' 'b'};
figure;
clear ULR SLR
for e = 1:numel(LR)
    ULR(e,:) = mean(LR{e},1);
    SLR(e,:) = std(LR{e},1,1)*size(LR{e},1).^-.5;
    errorbar(ULR(e,:),SLR(e,:),CL{e})
    hold all
end
errorbar(CONU,CONS,'k','LineWidth',3)
figure
hold on
for e = 1:numel(LR)
    plot(ULR(e,:),CL{e});
end
plot(CONU,'g','LineWidth',3)
