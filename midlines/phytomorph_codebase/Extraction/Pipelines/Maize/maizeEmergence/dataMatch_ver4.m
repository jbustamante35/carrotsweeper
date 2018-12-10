%% render mat files for movie measures
FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% extract
%O = cell(1,numel(FileList));
O = {};
w = [];
parfor e = 1:numel(FileList)
    %if ~isempty(strfind(FileList{e},'20170411_Camera2'))
        tic
        tmp = load(FileList{e});
        M = tmp.miniStack;
        [sm,mm] = squareMini(M,20);
        [O{e},~,w(e)] = measureMovie(sm,mm,0);
        toc
        toc*numel(FileList)/10/60
    %end
end
%% pour into container
for e = 1:numel(nO)
    if ~isempty(nO{e})
        O{e} = nO{e};
        w(e) = nw(e);
        e
    end
end
%% try to pick on
tic
tmp = load(mNM{midx(2)});
M = tmp.miniStack;
[sm,mm] = squareMini(M,20);
testMea = measureMovie(sm,mm,0);
toc*numel(FileList)/10/60
figure
plot(testMea)
%%
O{midx(3)} = testMea(1:277);
%% view with overlay
close all
sel = 70;
sel = 1209;
viewData(FileList{sel},[],1,0,1);
%% read hand measurements for emergence
%em = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/Emergence_hand_score_merged.csv');
ml = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/EmergenceAssay_Master_list.csv');
%% make labels for honey comb NEED
import java.util.Map.*;
import java.util.HashMap;
LT = HashMap();
trK = HashMap();

LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
HL = {'d' 'p'};
LABELS = {};
cnt = 1;
for e1 = 1:numel(HL)
    for e2 = 1:numel(LL)
        for e3 = 1:numel(NL)
            LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
            rasterP = num2str(cnt);
            cnt = cnt + 1;
            trK.put(rasterP,LABELS{end});
        end
    end
end
%% pull data from new measures
import java.util.Map.*;
import java.util.HashMap;
FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
LT = HashMap();
FNK = HashMap();
for e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    fidx = strfind(FileList{e},'_');
    gidx = strfind(FileList{e},'.');
    tmpK = FileList{e}((fidx(end)+1):(gidx(end)-1));
    newN = trK.get(tmpK);
    newDataKey{e} = strrep(nm,['_' tmpK],['_' newN]);
    fidx = strfind(newDataKey{e},'_');
    newDataKey{e}(fidx(end)) = '-';
    
    
    tmpD = O{e}';
    %{
    for k = 1:size(tmpD,2)
        tmpD(:,k) = imfilter(tmpD(:,k),fspecial('average',[11 1]),'replicate');
        %tmpD(:,k) = subBL(tmpD(:,k),50);
    end
    %}
    if ~isempty(strfind(lower(newDataKey{e}),'20170411'))
        %break
        %e
        %newDataKey{e}
    end
    if strcmp(lower(newDataKey{e}),'20170411')
        break
        %tmpD
    end
    if size(tmpD,2) == 5
        LT.put(lower(newDataKey{e}),tmpD);

        FNK.put(lower(newDataKey{e}),FileList{e});
    else
        e
    end
end
%% align data from mastersheet STACKY1 for hand scores
close all
iplantName = 27+1;
hcName = 20+1;
posName = 10+1;
genoName = 7+1;
FRAME_SCORE = 34;
OFFSET = {};
OFFSET_W = {};
import java.util.Map.*;
import java.util.HashMap;
GT = HashMap();
GT2 = HashMap();
MT = HashMap();
GV = HashMap();
SM2 = HashMap();
X = [];
X2 = [];
X3 = [];
X4 = [];
Y = [];
Y2 = [];
Z = [];
Z1 = [];
Z2 = [];
Z3 = [];
KJ = {};
mV = {};
mmmV = {};
WHATEVER = {};
PLOTGN = {};
FR = 200;
oIDX = [];
for e = 2:size(em,1)
    try
        if ~isnan(em{e,FRAME_SCORE})
            gidx = strfind(em{e,iplantName},filesep);
            key = [lower(em{e,iplantName}((gidx(end)+1):(end))) '-' lower(em{e,hcName}(1)) '-' lower(em{e,29}(1)) lower(num2str(em{e,30}))];
            key
            %e
            if strcmp(key,lower('20170411_Camera1-d-A1'))
               break
            end
            value = LT.get(key);

            if ~isempty(value)% && ~ischar(em{e,FRAME_SCORE})%&& isempty(strfind(key,'20170220_camera2-')) 
                
                % stack the popped value
                %Y2 = [Y2 em{e,31}];
                %em{e,FRAME_SCORE}
                if (em{e,FRAME_SCORE} > 218)
                    POP = 0;
                else
                    POP = 1;
                end
                Y2 = [Y2 POP];
                % stack the data for popped prediction
                Z = [Z value(1:FR,1)];
                Z1 = [Z1 value(1:FR,2)];
                Z2 = [Z2 value(1:FR,3)];
                Z3 = [Z3 value(1:FR,5)];
                PLOTGN{end+1} = genoTmp;
                % get the genotype name
                genoTmp = em{e,genoName};
                % store the genotype under the key
                GT2.put(key,genoTmp);

                WHATEVER{end+1} = key;

                mmmV{end+1} = FNK.get(key);
                oIDX = [oIDX e-1];
                
                OFFSET_W{end+1} = -(toTM((em{e,15}),em{e,16}) - toTM((em{e,24}),em{e,25}))/60;
                
                %waitforbuttonpress
                if ~ischar(em{e,FRAME_SCORE})
                if ~isnan(em{e,FRAME_SCORE}) && ~isinf(em{e,FRAME_SCORE}) && em{e,FRAME_SCORE} ~= 0 && em{e,FRAME_SCORE} <= 288% && all(value(1:FR,3) ~= 0)
                    tmp = zeros(size(value,1),1);
                    tmp(em{e,FRAME_SCORE}:(em{e,FRAME_SCORE}+1)) = 1;%1*max(value);
                    tmp = cumsum(tmp)';
                    %{
                    % pulse vs heavyside
                    fidx = find(tmp==1);
                    tmp = zeros(size(tmp))';
                    tmp(fidx(1)) = 1;
                    %}
                    Y = [Y tmp(1:FR)'];
                    X = [X value(1:FR,1)];
                    X2 = [X2 value(1:FR,2)];
                    X3 = [X3 value(1:FR,3)];
                    X4 = [X4 value(1:FR,5)];

                    % store the key in order
                    KJ{end+1} = key;
                    mV{end+1} = FNK.get(key);
                    GT.put(key,genoTmp);
                    
                    OFFSET{end+1} = -(toTM((em{e,15}),em{e,16}) - toTM((em{e,24}),em{e,25}))/60;
                end
                end


            else
                key
                e
                value;
            end
        end
    catch ME
        ME
        e
    end
end

%% masterlist STACK
iplantName = 28;
hcName = 21;
posName = 11;
genoName = 8;
import java.util.Map.*;
import java.util.HashMap;
GT = HashMap();
MT = HashMap();
GV = HashMap();
SM2 = HashMap();
Z = [];
Z1 = [];
Z2 = [];
Z3 = [];
WHOLE_GENO = {};
WHATEVER = {};
OFFSET_W = {};
RTYPE = {};
IDXW = [];
for e = 2:size(ml,1)
    %for f = 1:numel(d)
        %if ~isempty(strfind(ml{e,iplantName},d(f).experimentName))
        if ~isempty(ml{e,iplantName})
            arrayName = ml{e,posName};
            arrayName = ['; ' arrayName ';'];
            fidx = strfind(arrayName,';');
            INCELL = 0;
            for l = 1:(numel(fidx)-1)
                
                snip = arrayName((fidx(l)+1):(fidx(l+1)-1));
                sidx = strfind(snip,':');
                letterRange1 = snip(2);
                letterRange2 = snip(sidx(1)+1);
                numberRange1 = str2num(snip(3:(sidx(1)-1)));
                numberRange2 = str2num(snip((sidx(1)+2):end));
                gidx = strfind(ml{e,iplantName},filesep);
                cn = numberRange1:numberRange2;
                letterRange1 = find(strcmp(LL,letterRange1));
                letterRange2 = find(strcmp(LL,letterRange2));
                rn = letterRange1:letterRange2;
                %culN = ml{e,6};
                
                for r = 1:numel(rn)
                    for c = 1:numel(cn)
                        key = [ml{e,iplantName}((gidx(end)+1):end) '-' lower(ml{e,hcName}(1)) '-' LL{rn(r)} num2str(cn(c))];
                        %flag  = ~strcmp(ml{e,6},'15A-6197-02');
                        
                        key = lower(key);
                        value = LT.get(key);
                        if ~isempty(value)
                            deltaT = find(value);
                            genoTmp = ml{e,genoName};

                            if strcmp(genoTmp,'Ark 0812-87ne')
                                INCELL = INCELL + 1;
                            end
                                
                            WHOLE_GENO{end+1} = genoTmp;
                            WHATEVER{end+1} = key;
                            Z = [Z value(1:FR,1)];
                            Z1 = [Z1 value(1:FR,2)];
                            Z2 = [Z2 value(1:FR,3)];
                            Z3 = [Z3 value(1:FR,5)];
                            IDXW = [IDXW e];
                            
                            OFFSET_W{end+1} = -(toTM((ml{e,15}),ml{e,16}) - toTM((ml{e,24}),ml{e,25}))/60;
                            key;
                            RTYPE{end+1} = ml{e,5};
                            %{
                            posTmp = EM.get(key);
                            vec = GT.get(genoTmp);
                            vec = [vec value(1:250)];
                            if flag
                                GT.put(genoTmp,vec);
                            end


                            if isempty(deltaT)
                                deltaT = inf;
                            end
                            corrV = [posTmp,deltaT(1)];
                            if flag
                                MT.put(key,corrV);
                            end

                            vec = GV.get(genoTmp);
                            vec = [vec corrV'];
                            if flag
                                GV.put(genoTmp,vec);
                            end

                            vec = SM2.get(culN);
                            vec = [vec;deltaT(1)];
                            if flag
                                SM2.put(culN,vec);
                            end
                            %}
                        end
                    end
                end
                
            %end
            lower(ml{e,hcName}(1));
            end
            if strcmp(genoTmp,'Ark 0812-87ne')
                e
                INCELL
            end
        end
    %end
end

%% train for germination
toMG = {};
toMG{end+1} = Z1;%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});
toMG{end+1} = Z3;%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});
toMG{end+1} = bsxfun(@times,Z1,max(Z1,[],1).^-1);%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = bsxfun(@times,Z3,max(Z3,[],1).^-1);%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = cumsum(toMG{end});%bsxfun(@minus,X,mean(X,1));
%%
ncg = [3 3 3 3 3 3 3 3];
gMM = [];
for e = 1:numel(toMG)
    [wS wCG{e} wUG{e} wEG{e}] = PCA_FIT_FULL_T(toMG{e},ncg(e));
    gMM = [gMM;wCG{e}];
end

sY2 = zeros(2,size(Y2,2));
sY2(1,find(Y2==0)) = 1;
sY2(2,find(Y2==1)) = 1;
clear netG;
c = [];
for e = 1:20
    netG{e} = patternnet(10);
    netG{e} = train(netG{e},gMM,sY2);
germ = sim(netG{e},gMM);

[c(e),cm(:,:,e),ind,per] = confusion(sY2,germ);
cm(:,:,e)
end
[mis,midx] = min(c);
cm(:,:,midx)

germNet = netG{midx};
germ = sim(germNet,gMM);

%% score germ train
cnt =1;
TR = linspace(0,1,100);
for th = TR
    
    as = germ(1,:) > th;
    as = [as;~as];
    
    max(as,[],1);
    vecH = zeros(2,size(as,2));
    vecA = zeros(2,size(as,2));
    
    
    vecH(1,find(Y2==0)) = 1;
    vecH(2,find(Y2==1)) = 1;
    
    %vecA(1,find(as==0)) = 1;
    %vecA(2,find(as==1)) = 1;
    vecA = as;
    [c(cnt),cm(:,:,cnt),ind,per] = confusion(vecH,vecA);
    cm
    cnt = cnt + 1;
end

%%
[~,midx] = min(c);
cm(:,:,midx)
%%
thresholdV = TR(midx);
as = germ(1,:) > TR(midx);
    as = [as;~as];
vecH = zeros(2,size(as,2));
vecA = zeros(2,size(as,2));


vecH(1,find(Y2==0)) = 1;
vecH(2,find(Y2==1)) = 1;

%vecA(1,find(as==0)) = 1;
%vecA(2,find(as==1)) = 1;
vecA = as;
[c,cm,ind,per] = confusion(vecH,vecA);
sum(germ<.5);
%%
fidx = find((vecH(2,:) == 1) & (vecA(1,:) == 1));
fidx = find((vecH(1,:) == 1) & (vecA(2,:) == 1));
viewData(mmmV{fidx(5)},[],1,0,0);
% 1 3 4 8 10
% 1 5
WHATEVER{fidx(5)}
%% STACKY0 for germination frame
toM = {};
%toM{1} = X5;%bsxfun(@minus,X,mean(X,1));
%toM{2} = gradient(toM{1});%bsxfun(@minus,X,mean(X,1));
%toM{2} = cumsum(toM{1});%bsxfun(@minus,X,mean(X,1));
toM{end+1} = X2;%X2;%bsxfun(@minus,X,mean(X,1));
toM{end+1} = gradient(toM{end});
toM{end+1} = X4;%X2;%bsxfun(@minus,X,mean(X,1));
toM{end+1} = gradient(toM{end});
toM{end+1} = bsxfun(@times,X2,max(X2,[],1).^-1);%X2;%bsxfun(@minus,X,mean(X,1));
toM{end+1} = gradient(toM{end});%bsxfun(@minus,X,mean(X,1));
%toM{end+1} = X3;%bsxfun(@minus,X,mean(X,1));
%toM{end+1} = gradient(toM{end});%bsxfun(@minus,X,mean(X,1));
toM{end+1} = bsxfun(@times,X4,max(X4,[],1).^-1);%bsxfun(@minus,X,mean(X,1));
toM{end+1} = cumsum(toM{end});%bsxfun(@minus,X,mean(X,1));

%{
toM{1} = X4;%bsxfun(@minus,X,mean(X,1));
toM{2} = gradient(toM{1});%bsxfun(@minus,X,mean(X,1));
%}
%{
toM{1} = X3.*X2;%bsxfun(@minus,X,mean(X,1));
toM{2} = gradient(toM{1});%bsxfun(@minus,X,mean(X,1));
%}
%toM{2} = X2;%bsxfun(@minus,X,mean(X,1));
%toM{2} = bsxfun(@times,X,max(X,[],1).^-1);
%toM{3} = gradient(toM{1});%bsxfun(@minus,X,mean(X,1));
%toM{1} = bsxfun(@times,X,max(X,[],1).^-1);
%toM{2} = X2;%bsxfun(@minus,X2,mean(X2,1));
%toM{3} = X3;%bsxfun(@minus,X3,mean(X3,1));
%toM{4} = cumsum(toM{1},1);
%toM{5} = cumsum(toM{2},1);
%toM{6} = cumsum(toM{3},1);
%[~,toM{4}] = gradient(X);
%[~,toM{5}]= gradient(X2);
%[~,toM{6}]= gradient(X3);
%% STACKY2_v2
nf = 13;
nc = 3;
clear wS wC wU wE;
gM = [];
nc = [5 5 5 5 5 5 4];
nc = [3 3 3 3 3 3 3 3 4];
%nc = [2 2 2 2 2 2 4];
clear wC
for e = 1:numel(toM)
    WIN = im2col(toM{e},[nf 1],'sliding');
    [wS wC{e} wUF{e} wEF{e}] = PCA_FIT_FULL_T(WIN,nc(e));
    
    gM = [gM;wC{e}];
end
% make the Y var
WINY = im2col(Y,[nf 1],'sliding');
WINY = WINY((end-1)/2,:);
frt = col2im(WINY,[nf 1],size(X),'sliding');

%{
% make the global measures
[Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,lambda] = plsregress(gM',WINY',nc(end));
%lambda = mynLDA(ob,WINY,1,3);
wC{end+1} = (bsxfun(@minus,gM,mean(gM,1))'*lambda)';
%}
%{
cl = zeros(3,size(WINY,2));
for e = 1:size(cl,1)
    cl(e,WINY==(e-1)) = 1;
end
%}
cl = WINY;
%net = feedforwardnet([5 10 8]);
net = feedforwardnet([10 5]);
%net = patternnet(10);
%net.layers{end}.transferFcn = 'logsig';
net = train(net,gM,cl,'useParallel','yes');
wC{end+1} = sim(net,gM);
%{
%net = patternnet(20);
net.layers{end}.transferFcn = 'logsig';
net = train(net,gM,cl,'useParallel','yes');
tmp = sim(net,gM);
[~,wC{end+1}] = max(tmp,[],1);
%wC{end} = tmp;
%}


SIG = {};
for e = 1:numel(wC)
    SIG{e} = [];
    for t = 1:size(wC{e},1)
        tmp = col2im(wC{e}(t,:),[nf 1],size(X),'sliding');
        %tmp = imfilter(tmp,fspecial('average',[31 1]),'replicate');
        SIG{e} = cat(3,SIG{e},tmp);
    end
end
%WINY= cumsum(frt,1);
%WINY = WINY(:);
%{
for e = 1:(size(numel(wC)-1))
    for e = 1:size(wC{e},1)
        SIG{e}(:,:,e) = bsxfun(@times,SIG{e}(:,:,e),SIG{end});
    end
end
%}
%% LOOKUP PICK ON
for e = 1:numel(FileList)
    G(e) = ~isempty(strfind(FileList{e},'20170411_Camera2_8.mat'));
end
keyL = lower('20170411_camera2-d-a8');
GG = find(strcmp(WHATEVER,keyL));
germ(GG)
WHATEVER(GG)
FileList{G}
germ(GG)
K = find(frame(:,GG));
K(1)
%% pick on
viewData(FileList{find(G)},[],1,0,0);

%% SCORE ALL



toMG = {};
toMG{end+1} = Z1;%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});
toMG{end+1} = Z3;%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});
toMG{end+1} = bsxfun(@times,Z1,max(Z1,[],1).^-1);%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = bsxfun(@times,Z3,max(Z3,[],1).^-1);%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = cumsum(toMG{end});%bsxfun(@minus,X,mean(X,1));

toM = toMG;


[germ frame] = score(toMG,toM,wUG,wEG,wUF,wEF,germNet,net,nf,hmm,thresholdV);
%germ = ~germ;
%frameBK = frame;
frame(:,germ==0) = 0;
%{
%% align genotype SCORE2
for e = 1:numel(subJK)
    geno{e} = GT.get(subJK{e});
    
end
clear para
close all
%%
geno = PLOTGN;
%}
%% fit data
close all
%UQ = unique(geno);
UQ = unique(WHOLE_GENO);
para = {};
ALGO = [];
ALGOf = [];
X = 1:12*24*60;
subTYPE = {};
for u = 1:1%numel(UQ)
    %{
    fidx = strcmp(geno,UQ{u});
    subD = algoD(:,fidx);
    binA = mean(subD,2);
    subD = handD(:,fidx);
    binH = mean(subD,2);
    binH = cumsum(binH,1);
    %}
   
    
    fidx = strcmp(WHOLE_GENO,UQ{u});
    tmpO = cell2mat(OFFSET_W(find(fidx)));
    
    subKEY = WHATEVER(fidx);
    FIL = [];
    for e = 1:numel(subKEY)
        FIL(e) = logical(strcmp(subKEY{e}(1:16),'20170317_camera2'));
        FIL(e) = FIL(e) | logical(strcmp(subKEY{e}(1:18),'20170411_camera1-p'));
        FIL(e) = FIL(e) | logical(strcmp(subKEY{e}(1:18),'20170411_camera3-d'));
    end
    
    subD = frame(:,fidx);
    subD(:,logical(FIL)) = [];
    if strcmp(UQ{u},'Mo18W')
        %break
    end
    subDD = zeros(numel(X),size(subD,2));
    for tr = 1:size(subD,2)
        tmp = find(subD(:,tr)==1);
        if ~isempty(tmp)
            tmp = tmp(1)*30 + round(tmpO(tr));
            subDD(X(tmp):end,tr)= 1;
        end
    end
    NUM(u) = sum(subDD(end,:));
    TOT(u) = size(subDD,2);
    
    ALGO(:,u) = mean(subDD,2);
    %subD = handD(:,fidx);
    %binH = mean(subD,2);
    %HAND(:,u) = cumsum(binH,1);
    
    
    
    
    sig = ALGO(:,u)';
    tmpX = 1:size(ALGO,1);
    xlab = X;
    
    
    
    [J xval] = min(abs(sig - mean(sig)));
    [para{u}] = fminsearch(@(X)mySigmoid_ver0(xlab,X,sig),[sig(end) .5 xlab(xval)]);
    u
    [e,ALGOf(:,u)] = mySigmoid_ver0(xlab',para{u});
    fidx = find(fidx);
    subTYPE{end+1} = RTYPE{fidx(1)};
end
%% plot
close all
CL = {'r' 'g' 'b' 'c' 'k' 'y' 'r--' 'g--' 'b--' 'r' 'g' 'b' 'c' 'k' 'y' 'r--' 'g--' 'b--' 'r' 'g' 'b' 'c' 'k' 'y' 'r--' 'g--' 'b--' 'r' 'g' 'b' 'c' 'k' 'y' 'r--' 'g--' 'b--'};
%{
figure;
hold on
for u = 1:size(HAND,2)
    plot(HAND(:,u),CL{u});
    
end
%}


kpidx = find(strcmp(subTYPE,'Cotton'));
TOT = TOT(kpidx);
NUM = NUM(kpidx);
UQ = UQ(kpidx);
ALGO = ALGO(:,kpidx);
ALGOf = ALGOf(:,kpidx);
legend(UQ)
figure
hold on
for u = 1:size(ALGO,2)
    plot(ALGOf(:,u),CL{u});
    LEG{u} = [UQ{u} '-' num2str(ALGOf(end,u)) '-' num2str(NUM(u)) '-' num2str(TOT(u))];
end
legend(LEG)
r = gradient(ALGOf')';
[mr,fr] = max(r,[],1);
WID = [];
for e = 1:size(r,2)
    WID(e) = sum(r(:,e) > .70*mr(e));
    WID(e) = sum(ALGOf(:,e) < .9*ALGOf(end,e) & ALGOf(:,e) > .1*ALGOf(end,e) );
end
WID = WID / 60/24;
figure
hold on
for u = 1:size(ALGO,2)
    plot(r(:,u),CL{u});
    
end
legend(UQ);





H = bsxfun(@times,r,sum(r,1).^-1);
L = log(H);
L(isinf(L)) = 0;
E = -mean(L.*H,1);
figure;
plot(WID,E,'.')
figure;
plot(ALGOf(end,:),WID,'.')
figure;
plot(ALGOf(end,:),E,'.')
figure;
plot(E,WID,'.')

for e = 1:numel(UQ)
    pubSHEET{e,1} = LEG{e};
    pubSHEET{e,2} = ALGOf(end,e);
    pubSHEET{e,3} = mr(e);
    pubSHEET{e,4} = fr(e);
    pubSHEET{e,5} = WID(e);
end

cell2csv('/mnt/snapper/nate/mirror_images/maizeData/jgustin/pubSHEET_cotton.csv',pubSHEET);
%csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/fourMeasures.csv',[ALGOf(end,:));
%csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/bioD.csv',ALGOf);
%csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/biogD.csv',r);
%{
%% regress delta
close all
delta = abs((tr(:,1) - tr(:,2)));
[S vecs U E L ERR LAM] = PCA_FIT_FULL_T(vec,3);
netF = feedforwardnet([10 5 10]);

%netF = feedforwardnet([5 10 5]);
%net = feedforwardnet([20]);
%net = patternnet(15);
%netF.layers{end}.transferFcn = 'logsig';
netF = train(netF,vecs,delta','useParallel','yes');
fixV = sim(netF,vecs);
figure;
plot(delta',fixV,'.')
%}
%% 
[J,sidx] = sort(delta,'descend');
subJK(sidx(1:30))
%%
rm = fixV > 15;
tr(rm,:) = [];
%%
close all
plot(tr(:,1),tr(:,2),'.')
corr(tr(:,1),tr(:,2))
mean(tr(:,1)-tr(:,2))
std(tr(:,1)-tr(:,2))
%%
%{
%% STACKY2
nf = 21;
nc = 5;
WIN = im2col(X,[nf 1],'sliding');
[wS wC wU wE] = PCA_FIT_FULL_T(WIN,nc);

WIN2 = im2col(X2,[nf 1],'sliding');
[wS2 wC2 wU2 wE2] = PCA_FIT_FULL_T(WIN2,nc);

[~,gX] = gradient(X);
WIN3 = im2col(gX,[nf 1],'sliding');
[wS3 wC3 wU3 wE3] = PCA_FIT_FULL_T(WIN3,nc);

[~,gX2] = gradient(X2);
WIN4 = im2col(gX2,[nf 1],'sliding');
[wS4 wC4 wU4 wE4] = PCA_FIT_FULL_T(WIN4,nc);



WINY = im2col(Y,[nf 1],'sliding');
WINY = WINY((end-1)/2,:);
frt = col2im(WINY,[nf 1],size(X),'sliding');


ob = [wC;wC2;wC3;wC4]';
[Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,lambda] = plsregress(ob,WINY',nc);
%lambda = mynLDA(ob,WINY,1,3);
wC5 = (bsxfun(@minus,ob,mean(ob,1))*lambda)';
%wC4 = col2im(wC4',[nf 1],size(X),'sliding');
%{
tmpC = wC3;
wC3 = wC4;
wC4 = tmpC;
%}
SIG = [];
SIG2 = [];
SIG3 = [];
SIG4 = [];
SIG5 = [];
for e = 1:size(wC,1)
    tmp = col2im(wC(e,:),[nf 1],size(X),'sliding');
    tmp2 = col2im(wC2(e,:),[nf 1],size(X2),'sliding');
    tmp3 = col2im(wC3(e,:),[nf 1],size(gX),'sliding');
    tmp4 = col2im(wC4(e,:),[nf 1],size(gX),'sliding');
    tmp5 = col2im(wC5(e,:),[nf 1],size(gX),'sliding');
    %tmp = sum(tmp,1);
    SIG = cat(3,SIG,tmp);
    SIG2 = cat(3,SIG2,tmp2);
    SIG3 = cat(3,SIG3,tmp3);
    SIG4 = cat(3,SIG4,tmp4);
    SIG5 = cat(3,SIG5,tmp5);
end
%}
%% build
sig = mean(frt,2)';
xlab = 1:numel(sig);
[J xval] = min(abs(sig - mean(sig)));
[para] = fminsearch(@(X)mySigmoid_ver0((1:numel(sig)),X,sig),[sig(end) .5 xlab(xval)]); 
[e,sigF] = mySigmoid_ver0(xlab',para);
myT = sigmoidTransitionFunction(para,0);
myTi = sigmoidTransitionFunction(para,1);
%% basic massive
clear postU preU postCov preCov preDis postDis
preIdx = find(WINY==0);
tranIdx = find(WINY==1);
postIdx = find(WINY==2);

preGerm = hmm_node('preGerm');
tranGerm = hmm_node('tranGerm');
postGerm = hmm_node('postGerm');



for m = 1:numel(wC)
    
    preCov{1,m} = cov(wC{m}(:,preIdx)');
    preU{1,m} = mean(wC{m}(:,preIdx)');
    
    
    tranCov{2,m} = cov(wC{m}(:,tranIdx)');
    tranU{2,m} = mean(wC{m}(:,tranIdx)');
    
    postCov{3,m} = cov(wC{m}(:,postIdx)');
    postU{3,m} = mean(wC{m}(:,postIdx)');
    
    preDis{m} = myProb(preU{1,m},preCov{1,m});
    tranDis{m} = myProb(tranU{2,m},tranCov{2,m});
    postDis{m} = myProb(postU{3,m},postCov{3,m});
    %{
    preDis{m}.fitToGMM(wC{m}(:,preIdx)',2,1);
    postDis{m}.fitToGMM(wC{m}(:,postIdx)',2,1);
    %}
    %{
    preDis{m}.fitToGMM(wC{m}(:,preIdx)',2,m);
    tranDis{m}.fitToGMM(wC{m}(:,postIdx)',2,m);
    %}
    postDis{m}.fitToGMM(wC{m}(:,postIdx)',2,m);
    
    
    
    preGerm.attachDistribution(preDis{m},m);
    tranGerm.attachDistribution(tranDis{m},m);
    postGerm.attachDistribution(postDis{m},m);
end

%{

preCov2 = cov(wC2(:,preIdx)');
preU2 = mean(wC2(:,preIdx)');
postCov2 = cov(wC2(:,postIdx)');
postU2 = mean(wC2(:,postIdx)');

preCov3 = cov(wC3(:,preIdx)');
preU3 = mean(wC3(:,preIdx)');
postCov3 = cov(wC3(:,postIdx)');
postU3 = mean(wC3(:,postIdx)');

preCov4 = cov(wC4(:,preIdx)');
preU4 = mean(wC4(:,preIdx)');
postCov4 = cov(wC4(:,postIdx)');
postU4 = mean(wC4(:,postIdx)');

preCov5 = cov(wC5(:,preIdx)');
preU5 = mean(wC5(:,preIdx)');
postCov5 = cov(wC5(:,postIdx)');
postU5 = mean(wC5(:,postIdx)');
%}
preTopre = constantTransitionFunction(.5);
preTotran = constantTransitionFunction(.5);
%tranTotran = constantTransitionFunction(.5);
%t1 = heavisideTransitionFunction(20,@(x,y)lt(x,y));
tranTopost = constantTransitionFunction(1);
postTopost = constantTransitionFunction(1);


preGerm.attachNode(preGerm,preTopre);
preGerm.attachNode(tranGerm,preTotran);
tranGerm.attachNode(postGerm,tranTopost);
postGerm.attachNode(postGerm,postTopost);
%preGerm.attachNode(preGerm,myTi);
%preGerm.attachNode(postGerm,myT);
%postGerm.attachNode(postGerm,postTopost);

%{
preDis = myProb(preU,preCov);
postDis = myProb(postU,postCov);
preDis2 = myProb(preU2,preCov2);
postDis2 = myProb(postU2,postCov2);
preDis3 = myProb(preU3,preCov3);
postDis3 = myProb(postU3,postCov3);
preDis4 = myProb(preU4,preCov4);
postDis4 = myProb(postU4,postCov4);
preDis5 = myProb(preU5,preCov5);
postDis5 = myProb(postU5,postCov5);

%}

%{
preGerm.attachDistribution(preDis,1);
preGerm.attachDistribution(preDis2,2);
preGerm.attachDistribution(preDis3,3);
preGerm.attachDistribution(preDis4,4);
preGerm.attachDistribution(preDis5,5);

postGerm.attachDistribution(postDis,1);
postGerm.attachDistribution(postDis2,2);
postGerm.attachDistribution(postDis3,3);
postGerm.attachDistribution(postDis4,4);
postGerm.attachDistribution(postDis5,5);

%}




hmm = my_hmm();
hmm.addNode(preGerm);
hmm.addNode(tranGerm);
hmm.addNode(postGerm);
hmm.dn = ones(numel(wC),1);
%% train hmm
UD = {};
observation_labels = [];
for e = 1:numel(wC)
    observation_labels = [observation_labels;e*ones(size(wC{e},1),1)];
end
%observation_labels = [ones(size(wC,1),1);2*ones(size(wC2,1),1);3*ones(size(wC3,1),1);4*ones(size(wC4,1),1);5*ones(size(wC5,1),1)];
for e = 1:size(SIG{1},2)
    STIM = [];
    for s = 1:numel(SIG)
         tmpV = squeeze(SIG{s}(:,e,:))';
        STIM = [STIM;tmpV];
    end
   
    UD{e} = [frt(:,e)'+1;STIM];
end

hmm.update(UD,observation_labels);
%%
FIGVEC = [];
for e = 1:size(SIG{1},2)
STIM = [];
for s = 1:1%numel(SIG)
    tmpV = squeeze(SIG{s}(:,e,:))';
    STIM = [STIM;tmpV];
    
end
FIGVEC = [FIGVEC;STIM'];
e
end

%%
[fS fC fU fE fL fERR fLAM] = PCA_FIT_FULL(FIGVEC,2);
close all
figure;
for e = 1:size(SIG{1},2)
STIM = [];
for s = 1:1%numel(SIG)
    tmpV = squeeze(SIG{s}(:,e,:))';
    STIM = [STIM;tmpV];
    
end
    [ftC] = PCA_REPROJ(STIM',fE,fU);
    phase1 = find(frt(:,e)==0);
    phase2 = find(frt(:,e)==1);
    phase3 = find(frt(:,e)==2);
    plot(ftC(phase1,1),ftC(phase1,2),'r.','MarkerSize',1)
    hold on
    plot(ftC(phase2,1),ftC(phase2,2),'g.','MarkerSize',1)
    plot(ftC(phase3,1),ftC(phase3,2),'b.','MarkerSize',1)
    hold on 
    drawnow
end
%% apply HM
close all
tr = [];
subK = {};
IDX = [];
n = figure;
cor = figure;
mon = figure;
CV = [];
mea = [];
GT = [];
AT = [];
stdM = [];
vecGG = [];
mNM = {};
%observation_labels = [ones(size(wC,1),1);2*ones(size(wC2,1),1);3*ones(size(wC3,1),1);4*ones(size(wC4,1),1);5*ones(size(wC5,1),1)];
%observation_labels = [ones(size(wC,1),1);2*ones(size(wC2,1),1);3*ones(size(wC3,1),1);4*ones(size(wC4,1),1)];
%observation_labels = [ones(size(wC,1),1);2*ones(size(wC2,1),1);3*ones(size(wC3,1),1)];
%observation_labels = [ones(size(wC,1),1);2*ones(size(wC2,1),1)];
for e = 1:size(SIG{1},2)
    
    STIM = [];
    for s = 1:numel(SIG)
        tmpV = squeeze(SIG{s}(:,e,:))';
        STIM = [STIM;tmpV];
    end
    
    %{
    tmpV = squeeze(SIG(:,e,:));
    tmpV2 = squeeze(SIG2(:,e,:));
    tmpV3 = squeeze(SIG3(:,e,:));
    tmpV4 = squeeze(SIG4(:,e,:));
    tmpV5 = squeeze(SIG5(:,e,:));
    %}
    gidx = hmm.Viterbi(STIM,observation_labels,1);
    
    
    
    
    
    %close all
    fidx = find(gidx==2);
    fidx2 = find(frt(:,e)==2);
    
    if ~isempty(fidx) & ~isempty(fidx2)
        tr = [tr;[fidx2(1) fidx(1)]];
        subK{end+1} = KJ{e};
        IDX = [IDX ; e];
        vecGG = cat(3,vecGG,[toM{1}(:,e)]);
        mNM{end+1} = mV{e};
        GT = [GT frt(:,e)];
        AT = [AT gidx'==2];
        
        
        figure(n);
        os = bindVec(toM{1}((nf-1)/2:(end-(nf-1)/2),e));
        os2 = bindVec(toM{3}((nf-1)/2:(end-(nf-1)/2),e));
        %os3 = bindVec(toM{5}((nf-1)/2:(end-(nf-1)/2),e));
        MAG = max([max(os)]);
        plot(MAG*frt(:,e),'b');hold all;plot(MAG*(gidx-1),'k');%plot(os*5000 > .05,'g');
        plot(os,'b');
        plot(os2,'g');
        %plot(os3,'');
        hold off
        title(subK{end});
        
        %{
        figure(mon);
        testM = viewMon(mNM{end},fidx2(1),fidx(1));
        testM = imresize(testM,1);
        imshow(testM,[])
        title([num2str(tr(end,1)) '--' num2str(tr(end,2))]);
        hold off
        %waitforbuttonpress
        %}
        title(subK{end});
    hold off
    
    
    
    e
    CV = [CV corr(tr(:,1),tr(:,2))];
    delta = tr(:,1) - tr(:,2);
    mea = [mea mean(delta)];
    stdM = [stdM std(delta)];
    figure(cor);
    errorbar(mea,stdM);
    drawnow
    
    else
        %break
    end
   
    
    
    
    %plot(os2,'r');
    %plot(os3,'g');
    if e < 20
        %waitforbuttonpress
    end
    
    drawnow
    
    
corr(tr(:,1),tr(:,2))
    %waitforbuttonpress
end
figure;
plot(tr(:,1),tr(:,2),'.')
corr(tr(:,1),tr(:,2))
%%
csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/tr.csv',tr);
%%
kp = [];
for e = 1:numel(subK)
    
    if ~isempty(strfind(subK{e},'20170220_camera4'))
        kp(e) = 1;
    else
        kp(e) = 1;
        
    end
    kp(e) = 1;
end
kp = logical(kp);
figure;
plot(tr(kp,1),tr(kp,2),'.')
corr(tr(kp,1),tr(kp,2))
corr(tr(:,1),tr(:,2))
deltaA = abs(tr(:,1) - tr(:,2));
delta = tr(kp,1) - tr(kp,2);
mean(deltaA)
mean(delta)
deltaA = abs((tr(:,1)-mean(delta)) - tr(:,2));

mean(deltaA)
%% find one

%% sub intestigate
close all
kidx = find(kp);
delta = tr(kp,1) - tr(kp,2);
[J,midx] = sort(abs(delta),'descend');


sel = 200;
badKey = subK{kidx(midx(sel))};
sIDX = IDX(kidx(midx(sel)));
sig = squeeze(vecGG(:,:,kidx(midx(sel))));
tmp = [];
for e = 1:numel(toM)
    tmp = [tmp,bindVec(toM{e}(:,sIDX))];
end
sig = toM{1}(:,sIDX);
sigToSave = [toM{1}(:,sIDX),toM{3}(:,sIDX)];
csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/twoSigs.csv',sigToSave);
%sig = tmp;
figure;
plot(sig);
title(badKey);
hand = zeros(size(sig));
algo = zeros(size(sig));
hand(tr(kidx(midx(sel)),1)) = 1;
algo(tr(kidx(midx(sel)),2)) = 1;
hold on
plot(algo*max(sig(:)),'k');
plot(hand*max(sig(:)),'b');

deltaA = abs(tr(kidx(midx(4:end)),1) - tr(kidx(midx(4:end)),2));
delta = tr(kidx(midx(4:end)),1) - tr(kidx(midx(4:end)),2);
mean(deltaA)
mean(delta)

deltaA = abs((tr(kidx(midx(4:end)),1)-mean(delta)) - tr(kidx(midx(4:end)),2));

mean(deltaA)
%waitforbuttonpress
%close all
%% toRemove
rm = kidx(midx([4]));
for e = 1:numel(subK)
    if ~isempty(strfind(subK{e},'20170220_camera2-'));
        rm = [rm e];
    end
end
trP = (tr*30 + [cell2mat(OFFSET(IDX))' cell2mat(OFFSET(IDX))'])/60;
trP(rm,:) = [];
corr(trP(:,1),trP(:,2))
close all
plot(trP(:,1),trP(:,2),'.')
delta = (trP(:,1)-trP(:,2));

mean(abs(delta))
mean(delta)

std(delta);
std([-2 2]);


csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/trFinal.csv',trP);
%%

%% view data
viewData(mNM{kidx(midx(sel))},sig',1,0,0);
%% JUNK

%{
%% look at NN vs real
close all
figure;
tr = [];
disp = 0;
earlyT = 40;
sc = [];
vec = [];
subJK = {};
algoD = [];
handD = [];
mNM = {};
vecGG = [];
for e = 1:size(frt,2)
    sig = SIG{end}(:,e);
    
    
    sig = sig / max(sig);
    sig = imfilter(sig,fspecial('average',[5 1]),'replicate');
    sig = sig / max(sig);
    osig = sig;
    pidx = find(sig == imdilate(sig,strel('disk',5)));
    sig(1:earlyT) = 0;
    ridx = sig(pidx) < .2;
    pidx(ridx) = [];
    
    plot(sig,'k');
    %{
    %fidx1 = find(SIG{end}(:,e) > .6);
    [~,fidx1] = max(SIG{end}(:,e));
    
    %}
    
    fidx = find(frt(:,e) == 1);
    tmp = zeros(size(frt(:,e)));
    tmp(pidx(1)) = 1;
    
    
    
    
    if ~isempty(fidx) && ~isempty(pidx)
        tr = [tr;[fidx(1) pidx(1)]];
        al = zeros(size(frt(:,1),1),1);
        al(pidx(1):end) = 1;
        %sc = [sc;sum(SIG{end}(:,e),1)];
        sc = [sc;[sum(osig,1) std(osig,1,1) sum(SIG{end}(:,e) < .003) sum(sig)]];
        vec = [vec sort(SIG{end}(:,e))];
        subJK{end+1} = KJ{e};
        algoD = [algoD al];
        handD = [handD frt(:,e)];
        mNM{end+1} = mV{e};
        vecGG = cat(3,vecGG,[toM{1}(:,e) toM{2}(:,e) toM{3}(:,e)]);
    end
    if disp
        hold on
        plot(frt(:,e)*5,'r')
        plot(tmp*5,'k');

        plot(toM{1}(:,e),'r');
        plot(toM{2}(:,e),'g');
        plot(toM{3}(:,e),'b');
        %{
        plot(toM{4}(:,e),'m');
        plot(toM{5}(:,e),'g');
        plot(toM{6}(:,e),'c');
        %}
        title(subJK{end});
        hold off
        %waitforbuttonpress
        drawnow
    end
    e
end
%sc = sum(SIG{end},1);
delta = abs(tr(:,1) - tr(:,2));
figure;
plot(sc(:,1),delta,'.')
figure;
plot(sc(:,2),delta,'.')
figure;
plot(prod(sc,2),delta);
figure;
plot(delta,sc(:,end))
figure
plot(tr(:,1),tr(:,2),'.')
corr(tr(:,1),tr(:,2))
mean(tr(:,1)-tr(:,2))
std(tr(:,1)-tr(:,2))
csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/tr.csv',tr);
%}
%{
%% TRACK
close all
%[~,TR1] = max(frame,[],1);
%[~,TR2] = max(handD,[],1);
deltaA = abs(tr(:,1) - tr(:,2));
[J,midx] = sort(deltaA,'descend');
viewData(mNM{midx(1)},vecGG(:,:,midx(1))');
%}

