%% 0) load shape data
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
S = schemas(conn);
T = tables(conn);
P = columns(conn);
sid = 'kernel_3d';
for e = 1:size(P,1)
    if strcmp(P{e},sid)
        cNames = P{e,2};
        sidNUM = e;
    end
end
%% 0.1 construct query
q = ['select * from kernel_3d '...    
    'join kernels on kernel_3d.kernel_id = kernels.id ' ...
    'join kernel_plates on kernels.plate_id = kernel_plates.id '...
    'join population_lines on kernel_plates.population_line_id = population_lines.id'];
results = fetch(conn,q);
%% 0.2 put kernel features into hashmap
kernelVec = HashMap();
ROWnum = 6;
COLnum = 8;
ROWvec = {'A','B','C','D','E','F'};
COLvec = {'1','2','3','4','5','6','7','8'};
for e = 1:size(results,1)
    
    try
        NUMpos = results{e,45}-1;
        
        if NUMpos >= ROWnum*COLnum
            NUMpos = NUMpos-ROWnum*COLnum;
            LETTER = 'B';
        else
            LETTER = 'A';
        end
        
        rowN = mod(NUMpos,ROWnum) + 1;
        rowV = ['_' ROWvec{rowN}];

        colN = floor(NUMpos/ROWnum) + 1;
        colV = ['_' COLvec{colN}];
        
        
    catch ME

        LETTER = '*';
        
    end
    
    
    wellN = [rowV colV];
    
    
    plateN = results{e,50};
    key1 = [plateN '*' wellN];
    key2 = [plateN LETTER '*' wellN];    
    kernelVec.put(key1,results(e,1:42));
    kernelVec.put(key2,results(e,1:42));
end
%% 1) load csv masterlist
[masterList, result]= readtext('/mnt/spaldingdata/florida/files/nir_master_list.csv');
%% 2) load oil strach protien
[pC]= readtext('/mnt/spaldingdata/nate/COEFF.csv');
pC = cell2mat(pC);
pC(:,end) = [];
for i = 1:size(pC,1)
    npC(i,:) = pC(i,:)/norm(pC(i,:));
end
%% 3) make platename,genotype key value store from masterlist
import java.util.Map.*;
import java.util.HashMap;
K1 = 7;
V1 = 11;
P2G = HashMap();
for e = 1:size(masterList,1)
    key = masterList{e,K1};
    value = masterList{e,V1};    
    if ~isempty(value)
        P2G.put(key,value);
    end
end
%% 4) dig for tip angle data
FilePath = '/mnt/scratch5/takeshi_dev_run5/csvSPOOL/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% 5) load the plate name data translate list
[d,r] = readtext('/mnt/spaldingdata/nate/plateList.csv');
% parse the plate name data
for i = 1:size(d,1)
    d{i,1} = d{i,1}(1:end-2);
    if isempty(d{i,2})
        key = 'na';
    else
        key = d{i,2};
    end
    value = d{i,1};
    MAP.(key) = value;
end
%% 6) open tip angles
cnt = 1;
DT = [];
for i = 1:numel(FileList)
    [p n ext] = fileparts(FileList{i});
    % snip to folder
    fidx = strfind(n,'--');
    newName = n(fidx(end) + 2:end);
    fidx = strfind(newName,'_');
    if isfield(MAP,newName)
        plateName = MAP.(newName);
        % factor
        wellNames = newName(fidx(2)+1:end);
        D = csvread(FileList{i});
        wellNames = ['_' wellNames '_'];
        fidx = strfind(wellNames,'_');

        DT = [DT D'];
        for i = 1:numel(fidx)-1
            WNT{cnt} = wellNames(fidx(i)+1:fidx(i+1)-1);
            PLT{cnt} = plateName;
            cnt = cnt + 1;
        end
    else
        plateName = newName(1:fidx(2)-1);
    end
end
%% 7) clean tip angles
% clean via derivative
ridx = find(any(abs(diff(DT,1,1)) > 20/180*pi,1));
DT(:,ridx) = [];
PLT(ridx) = [];
WNT(ridx) = [];
% clean via total bend
ridx = find(any(abs(DT(1,:) - DT(end,:)) < 25/180*pi,1));
DT(:,ridx) = [];
PLT(ridx) = [];
WNT(ridx) = [];
%% 8) run import init from mongoDB
import phytoG.locked.Bpersist.Bos.implementations.*
import com.mongodb.*;
import java.util.Map.*;
import java.util.HashMap;
import phytoG.locked.BdataObjects.fileSystem.implementations.imageList;
import phytoG.locked.BdataObjects.BbioObjects.maize.spectraData;
oStore = OStore_mdb();
oStore.accessResource();
oStore.setCollection('specData');
%% 9) import spec and pull kernel shape data from hashmap
f.specData = [];
f.tipAngle = [];
f.genoType = {};
for e = 1:numel(WNT)
    qMap = HashMap();
    pn = PLT{e};
    qMap.put('_k_m._pnode._k_m._plateName',pn);
    wn = ['_' WNT{e}(1) '_' WNT{e}(2)];
    qMap.put('_k_m._pnode._k_m._wellName',wn);
    cursor = oStore.search(qMap);
    itr = cursor.iterator();
    if itr.size()
        n = itr.next();
        n = spectraData(n);
        spec = n.getSpectrum();

        specData = [];
        for i = 1:spec.size()
            specData(i,1) = str2num(spec.get(i-1));
        end
        e
        f.specData = [f.specData specData];
        f.tipAngle = [f.tipAngle DT(:,e)];
        f.genoType{end+1} = P2G.get(pn);
       %{
        sKey = [pn '*' wn];
        kVec = kernelVec.get(sKey);
        if isempty(kVec)
            sKey = [pn 'A*' wn];
            kVec = kernelVec.get(sKey);
        end
        %}
        
    end
end
f.specData(1:3,:) = [];
%% 10) hold out
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
tsS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',6);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    tsS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,tsS,'r');    
    drawnow
    
end
%% 11) create tip angle from spec data
close all
disp = 1;
toFitY = f.tipAngle';
toFitX = f.specData';
for L = 1:3
    for num = 6

        [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,5);    

        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);


        %{
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u2 = mean(toFitY(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s2 = std(toFitY(fidx,:),1,1);

            dif(u) = mean(u1*u2');
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            %}
        end
        DIF(num) = mean(dif);
    end
end
%% 11.2) create tip angle from spec data with hold out
close all
disp = 1;
toFitY = f.tipAngle';
toFitX = f.specData';
perDraw = .8;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
for L = 1:3

    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,15);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,5);    

    %%%%%%%%%%%%%%%
    % perform corr
    [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
    [mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
    % predict with test
    k = sC(Test==1,:)*A*inv(B);
    mk = sC(Test==1,:)*mA*inv(mB);


    k = PCA_BKPROJ(k,tE,tU);
    mk = PCA_BKPROJ(mk,tE,tU);
    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)

        fidx = find(strcmp([f.genoType(Test==1)],UQ{u}));
        wfidx = find(strcmp([f.genoType],UQ{u}));


        u1 = mean(k(fidx,:));
        mu1 = mean(mk(fidx,:));
        u2 = mean(toFitY(wfidx,:),1);


        s1 = std(k(fidx,:),1,1);
        ms1 = std(mk(fidx,:),1,1);
        s2 = std(toFitY(wfidx,:),1,1);


        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');
                hold on
                %errorbar(mu1*-180/pi,ms1*-180/pi,'b');
                errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                pause(1);
                hold off
            catch

            end
        end
   
    end
end
%% 11.3) create hold out via genotype
close all
disp = 1;
toFitY = f.tipAngle';
toFitX = f.specData';
UQ = unique([f.genoType]);
num = 15;
for L = 1:3
    for u = 1:numel(UQ)

        Test = strcmp([f.genoType],UQ{u});
        Train = ~ Test;


        [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,3);    

        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC(Train,:),tC(Train,:));
        % predict with test
        k = sC(find(Test==1),:)*A*inv(B);


        k = PCA_BKPROJ(k,tE,tU);

        u1 = mean(k);
        u2 = mean(toFitY(Test==1,:),1);

        s1 = std(k,1,1);
        s2 = std(toFitY(Test==1,:),1,1);

        
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                pause(1);
                hold off
            catch

            end
        end
    end
end
%% 11.7) create prediction plot(s) of PCA values for tip angle
% build holdout index
perDraw = .01;
toFitY = f.tipAngle';
toFitX = f.specData';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

num = 11;
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,5);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% 11.8) create oil, starch, protien from tip angles via CCA
% build holdout index
f.ospw = pC*f.specData;
num = 11;
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,num);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
for i = 1:size(f.ospw,1);
    figure;    
    plot(f.ospw(i,:),k(:,i),'.')
    title(TITLE{i})
    axis equal
end
%% 11.8) create oil, starch, protien from tip angles via CCA by genotype
% build holdout index
f.ospw = pC*f.specData;
num = 11;
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,num);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
UQ = unique([f.genoType]);
for i = 1:size(f.ospw,1);    
    figure;
    hold on
    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        realV = mean(f.ospw(i,fidx));
        preV = mean(k(fidx,i));        
        plot(realV,preV,'.');        
    end
    title(TITLE{i})
    axis equal
end
%% 11.8) create oil, starch, protien from tip angles via pls
% build holdout index
perDraw = .3;
f.ospw = pC*f.specData;
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
num = 11;
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,5);    

% pls regression
[XL,YL,XS,YS,BETA] = plsregress(tC,sC,5);
% predict
k = [ones(size(tC,1),1) tC]*BETA;


k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
for i = 1:size(f.ospw,1);
    figure;        
    %plot(f.ospw(i,Test==1),k(:,i),'.')
    plot(f.ospw(i,:),k(:,i),'.')
    title(TITLE{i})
    axis equal
end
%% 11.9) create oil, starch, protien from tip angles directly
% build holdout index
f.ospw = pC*f.specData;
toFitX = f.ospw';
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,4);


%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);

TITLE = {'oil','starch','protien','weight'};
UQ = unique([f.genoType]);
for i = 1:size(f.ospw,1);        
    
    xP = f.ospw(i,:)';
    yP = k(:,i);
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    %plot(tC(plotSet==1,i),k(:,i),'.')
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);

    
    title([TITLE{i} '----' num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% test myCCA
toFitY = f.tipAngle';
toFitX = f.specData';
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,5);
[A,B,r,U,V,stats] = canoncorr(sC,tC);
[wx,wy] = myCCA(sC,tC,3);
%%

close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
tsS = [];
fitX = f.specData';
fitY = f.tipAngle';
cnt = 1;
for num = 3:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        
        
        [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(fitX,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(fitY,3);
        
        
    %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(Train==1,yi),tmpX(Train==1,:));
            k = [k tmpX*b];
        end
        %}
        


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
        % predict
        %predictY = sC*A*inv(B);
        predictY = sC*A*B';



        %{
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}
        
        predictY = PCA_BKPROJ(predictY,tE,tU);
        
        
        trTOP = sum(predictY(Train==1,:).*fitY(Train==1,:),2);
        trBOTX = sum(predictY(Train==1,:).*predictY(Train==1,:),2);
        trBOTY = sum(fitY(Train==1,:).*fitY(Train==1,:),2);
        trCORR = trTOP.*trBOTX.^-.5.*trBOTY.^-.5;
        
        tsTOP = sum(predictY(Test==1,:).*fitY(Test==1,:),2);
        tsBOTX = sum(predictY(Test==1,:).*predictY(Test==1,:),2);
        tsBOTY = sum(fitY(Test==1,:).*fitY(Test==1,:),2);
        tsCORR = tsTOP.*tsBOTX.^-.5.*tsBOTY.^-.5;
        
        
        trainR(loop) = mean(trCORR);
        testR(loop) = mean(tsCORR);        
        
    end
    
    trU(cnt) = mean(trainR);
    trS(cnt) = std(trainR);
    
    tsU(cnt) = mean(testR);
    tsS(cnt) = std(testR);
    
    cnt = cnt + 1;
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,tsS,'r');    
    drawnow
    
end
%% try correlation via genotype
close all
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',5);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);


UQ = unique([f.genoType]);
tipFig = figure;
specFig = figure;
ccFig = figure;
ccFigT = figure;
ccFigS = figure;
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'}
for u = 1:numel(UQ)
    
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    
    %[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData(:,fidx)',5);
    %[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle(:,fidx)',3);
    
    
    % get the tip angle
    tmpTip = f.tipAngle(:,fidx)';
    % get the spec data
    tmpSpec = f.specData(:,fidx)';
    % get the means
    tU = mean(tmpTip,1);
    sU = mean(tmpSpec,1);
    % get the std
    tS = std(tmpTip,1,1)*size(tmpTip,1)^-.5;
    sS = std(tmpSpec,1,1)*size(tmpSpec,1)^-.5;
    
    if numel(fidx) > 1
        [A,B,r,U,V,stats] = canoncorr(sC(fidx,:),tC(fidx,:));
        %[A,B,r,U,V,stats] = canoncorr(sC,tC);
        
        
        
        
        figure(ccFigS);
        sData = PCA_BKPROJ(A(:,1)',sE,sU);
        plot(sData);
        hold all;
        
        figure(ccFigT);
        tData = PCA_BKPROJ(B(:,1)',tE,tU);
        plot(tData);
        hold all;
        
        
        LEG{u} = [UQ{u} '--' num2str(r(1)) '--' num2str(numel(fidx)) '--' num2str(stats.p(1))];
        figure(ccFig);
        scatter(U(:,1),V(:,1),CL{u});axis equal
        hold all
    else
        figure(ccFigT);
        scatter(0,0,CL{u});
        LEG{u} = [UQ{u} '---NA'];
    end
    
    
    figure(tipFig);
    errorbar(tU,tS);
    hold all
    
    figure(specFig);
    errorbar(sU,sS);
    hold all
    
    
end
figure(tipFig);
legend(UQ);
figure(specFig);
legend(UQ);
figure(ccFig);
legend(LEG);
figure(ccFigS);
legend(UQ);
figure(ccFigT);
legend(UQ);
%% try genotype cluster
figure;
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',5);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

UQ = unique([f.genoType]);
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
utC = [];
usC = [];
for u = 1:numel(UQ)
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    if numel(fidx) > 2
        usC = [usC;mean(sC(fidx,:))];
        utC = [utC;mean(tC(fidx,:))];
    end
end


[A,B,r,U,V,stats] = canoncorr(usC,utC);
sData = PCA_BKPROJ(A(:,1)',sE,sU);
tData = PCA_BKPROJ(B(:,1)',tE,tU);
figure;plot(sData);
scatter(U(:,1),V(:,1));axis equal;
axis([-3 3 -3 3]);
%%
f.ospw = pC*f.specData;
close all
disp = 1;
toFitY = f.tipAngle';
toFitX = f.specData';
%toFitX = f.ospw';
for num = 15
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,5);    
    
    
    [A,B,r,U,V,stats] = canoncorr(sC,tC);
    k = sC*A*inv(B);
    k = PCA_BKPROJ(k,tE,tU);



    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        u1 = mean(k(fidx,:));
        u2 = mean(toFitY(fidx,:),1);

        
        s1 = std(k(fidx,:),1,1);
        s2 = std(toFitY(fidx,:),1,1);
        
        dif(u) = mean(u1*u2');
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                %axis([0 61 -30 90])
                pause(1);
                hold off
            catch

            end
        end
        %}
    end
    DIF(num) = mean(dif);
end
sK = PCA_BKPROJ(A',sE,sU);
%{
for e = 1:size(k,1)
    plot(k(e,:),'k');
    hold on
    plot(f.tipAngle(:,e),'r');
    hold off
    drawnow
    pause(.4)
end
%}%% oil starch protien

f.ospw = pC*f.specData;

figure;
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

UQ = unique([f.genoType]);
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
utC = [];
usC = [];
for u = 1:numel(UQ)
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    if numel(fidx) > 2
        usC = [usC;mean(f.ospw(fidx,:))];
        utC = [utC;mean(tC(fidx,:))];
    end
end


[A,B,r,U,V,st
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
tsS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    tsS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,tsS,'r');    
    drawnow
    
endats] = canoncorr(usC,utC);
sData = PCA_BKPROJ(A(:,1)',sE,sU);
tData = PCA_BKPROJ(B(:,1)',tE,tU);
figure;plot(sData);

scatter(U(:,1),V(:,1));axis equal;
axis([-3 3 -3 3]);
%% please try
f.ospw = pC*f.specData;
close all
disp = 1;
toFitY = f.tipAngle';
toFitX = f.specData';
%toFitX = f.ospw';
for num = 15
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(toFitX,num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(toFitY,5);    
    
    
    [A,B,r,U,V,stats] = canoncorr(sC,tC);
    k = sC*A*inv(B);
    k = PCA_BKPROJ(k,tE,tU);



    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        u1 = mean(k(fidx,:));
        u2 = mean(toFitY(fidx,:),1);

        
        s1 = std(k(fidx,:),1,1);
        s2 = std(toFitY(fidx,:),1,1);
        
        dif(u) = mean(u1*u2');
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                %axis([0 61 -30 90])
                pause(1);
                hold off
            catch

            end
        end
        %}
    end
    DIF(num) = mean(dif);
end
sK = PCA_BKPROJ(A',sE,sU);
%{
for e = 1:size(k,1)
    plot(k(e,:),'k');
    hold on
    plot(f.tipAngle(:,e),'r');
    hold off
    drawnow
    pause(.4)
end
%}
%% size matters
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',size(f.specData',2));
newData = (diag(diag(sLAM).^-.5)*sC')';
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(newData,5);
%%
ng = 1;
kidx = kmeans(f.tipAngle',ng);
kidx = kmeans(f.specData',ng);
%% init first tip angle to zero
f.tipAngle2 = bsxfun(@minus,f.tipAngle,f.tipAngle(1,:));
%f.tipAngle2 = diff(f.tipAngle,1,1);
%% 
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',8);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',8);
T = [sC tC];
ng = 1;
kidx = kmeans(T,ng);
%% try upbulk 2 term
for i = 1:size(sC,2)
    for j = i:size(sC,2)
        sC = [sC sC(:,i).*sC(:,j)];
    end
end
sE%% explore nonlinear
plot3(sC(:,1),sC(:,2),tC(:,2),'.')
%% pls regression
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',11);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);
[XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
Y0 = [ones(size(sC,1),1) sC]*BETA;
for e = 1:size(Y0,1)
    plot(Y0(e,:),'r');
    sData = PCA_BKPROJ(A(:,1)',sE,sU);
    hold on
    plot(f.specData(:,e),'b');
    hold off
    pause(.3)
end
%% clustered cannon corr --> 
CL = {'r.','b.','g.','k.','m.'};
figure;
for num = 10
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);
    for e = 1:ng
        [A,B,r,U,V,stats] = canoncorr(sC(kidx==e,:),tC(kidx==e,:));
        sData = PCA_BKPROJ(A(:,1)',sE,sU);
        tData = PCA_BKPROJ(B(:,1)',tE,tU);
        plot(U(:,1),V(:,1),CL{e});axis equal
        hold on;
        R(num) = r(1);
    end
    hold off;
end

axis([-4 4 -4 4])
figure;plot(sData);
figure;plot(tData);
figure;plot(-f.tipAngle);
figure;plot(f.specData);
%%
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
tsS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    tsS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,tsS,'r');    
    drawnow
    
end


%%
pth = '/mnt/spaldingdata/Bessie/forPoster/';
csvwrite([pth 'tipAngle.csv'],f.tipAngle);
%%
csvwrite([pth 'corr.csv'],[U(:,1),V(:,1)]);
%% view via eye
OP = tData;
mea = OP*f.tipAngle;
[J sidx] = sort(mea);

hold on;
for e = 1:numel(sidx)
    plot(tS','b');
    hold on
    plot(tS(sidx(e),:),'r');
    hold off
    drawnow
end
%% predict spec
A = sC'/tC';
sim = A*tC';
M = PCA_BKPROJ(sim',sE,sU);
for e = 1:size(M,1)
    plot(M(e,:),'r');
    hold on
    plot(f.specData(:,e),'b');
    hold off
    pause(.3)
end
%% predict tip angle
A = tC'/sC';
sim = A*sC';
M = PCA_BKPROJ(sim',tE,tU);
for e = 1:size(M,1)
    plot(M(e,:),'r');
    hold on
    plot(f.tipAngle(:,e),'b');
    hold off
    pause(.3)
end
    %%
    cnt = 1;
    while (itr.hasNext())
        n = itr.next();
        n = spectraData(n);
        spec = n.getSpectrum();

        for i = 1:spec.size()
            specData(i,1) = str2num(spec.get(i-1));
        end
        DS(:,cnt) = specData;
        WNS{cnt} = n.getPlateName();
        PNS{cnt} = n.getWellName();
        cnt = cnt + 1
    end
end