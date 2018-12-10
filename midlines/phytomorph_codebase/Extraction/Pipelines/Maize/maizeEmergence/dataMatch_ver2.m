%% signal processing from hand measurements GEN d
FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
clear d;
cnt =1;
for e = 1:numel(FileList)
    fileName = FileList{e};
    fidx = strfind(fileName,'_');
    
    [pth,fn,ext] = fileparts(FileList{e});
    
    if ~isempty(strfind(fn,'Red'))
        tmp = readtext(FileList{e});
        d(cnt).redSig = cell2mat(tmp);
        
        tmpFileName = [pth filesep strrep(fn,'percentRed','stdColor') '.csv'];
        tmp = readtext(tmpFileName);
        d(cnt).ampSig = cell2mat(tmp);
        fn
        fidx = strfind(fn,'_');
        d(cnt).experimentName = fn(1:(fidx(end)-1))
        
        cnt = cnt +1;
    end
    
    
end
%%
FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
for e = 1:numel(FileList)
    tic
    tmp = load(FileList{e});
    M = tmp.miniStack;
    
    % HERE NEW JULY 13 2017
    [sm,mm] = squareMini(M,20);
    [O{e}] = measureMovie(sm,mm,0);
    toc*numel(FileList)/10/60
end



%% view with overlay
close all
sel = 70;
viewData(FileList{sel},O{sel});
%% read hand measuremetns for emergence
em = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/Emergence_hand_score_merged.csv');

%% align data from mastersheet STACKY1
close all
iplantName = 27+1;
hcName = 20+1;
posName = 10+1;
genoName = 7+1;
FRAME_SCORE = 34;
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
Y = [];
Y2 = [];
Z = [];
Z1 = [];
Z2 = [];
KJ = {};
mV = {};

for e = 2:size(em,1)
    if ~isnan(em{e,FRAME_SCORE})
        gidx = strfind(em{e,iplantName},filesep);
        key = [lower(em{e,iplantName}((gidx(end)+1):(end))) '-' lower(em{e,hcName}(1)) '-' lower(em{e,29}(1)) lower(num2str(em{e,30}))];
        if strcmp(key,'20170220_camera1-d-a7')
         %   break
        end
        value = LT.get(key);
        if ~isempty(value)
            
            
            Y2 = [Y2 em{e,31}];
            
            Z = [Z value(1:273,1)];
            Z1 = [Z1 value(1:273,2)];
            Z2 = [Z2 value(1:273,3)];
            
            genoTmp = em{e,genoName};
            GT2.put(key,genoTmp);
            
            
            
            tmp = zeros(size(value,1),1);
            
            %waitforbuttonpress
            if ~isinf(em{e,FRAME_SCORE}) && em{e,FRAME_SCORE} ~= 0 && em{e,FRAME_SCORE} <= 288
                
                tmp(em{e,FRAME_SCORE}:end) = 1;%1*max(value);
            
                % pulse vs heavyside
                fidx = find(tmp==1);
                tmp = zeros(size(tmp))';
                tmp(fidx(1)) = 1;
                
                Y = [Y tmp(1:273)'];
                %break
                
                %Y2 = [Y2 em{e,31}];
                X = [X value(1:273,1)];
                X2 = [X2 value(1:273,2)];
                X3 = [X3 value(1:273,3)];
                
                
                KJ{end+1} = key;
                mV{end+1} = FNK.get(key);
                
                
                GT.put(key,genoTmp);
                %{
            plot(value);
            hold on
            plot(tmp);
            hold off
            title(key);
                %}
            %axis([0 270 0 .0001]);
            %waitforbuttonpress
            end
            
            
            %Y2 = [Y2 em{e,FRAME_SCORE}];
            %drawnow
            %waitforbuttonpress
        else
            key
        end
    end
end
%% train for germ

toMG = {};
%oM{1} = X2;

toMG{1} = cumsum(Z,1);
toMG{2} = cumsum(Z1,1)
toMG{3} = cumsum(Z2,1);
ncg = [3 3 3];
gMM = [];
for e = 1:numel(toMG)
     %WIN = im2col(toMG{e},[nf 1],'sliding');
    [wS wCG{e} wUG{e} wEG{e}] = PCA_FIT_FULL_T(toMG{e},ncg(e));
    gMM = [gMM;wCG{e}];
end

sY2 = zeros(2,size(Y2,2));
sY2(1,find(Y2==0)) = 1;
sY2(2,find(Y2==1)) = 1;
clear netG;
for e = 1:40
    netG{e} = patternnet(10);
    netG{e} = train(netG{e},gMM,sY2);
germ = sim(netG{e},gMM);

[c(e),cm(:,:,e),ind,per] = confusion(sY2,germ);
cm(:,:,e)
end

germ = sim(netG{43},gMM);

%%
cnt =1;
TR = linspace(0,1,100);
for th = TR
    
    as = germ < th;
    vecH = zeros(2,size(as,2));
    vecA = zeros(2,size(as,2));
    
    
    vecH(1,find(Y2==0)) = 1;
    vecH(2,find(Y2==1)) = 1;
    
    vecA(1,find(as==0)) = 1;
    vecA(2,find(as==1)) = 1;
    
    [c(cnt),cm,ind,per] = confusion(vecH,vecA);
    cm
    cnt = cnt + 1;
end
%%
[~,midx] = min(c);
as = germ < TR(midx)
vecH = zeros(2,size(as,2));
vecA = zeros(2,size(as,2));


vecH(1,find(Y2==0)) = 1;
vecH(2,find(Y2==1)) = 1;

vecA(1,find(as==0)) = 1;
vecA(2,find(as==1)) = 1;

[c,cm,ind,per] = confusion(vecH,vecA);
sum(germ<.5)

%% NOPE
nf = 11;
WIN = im2col(X,[nf 1],'sliding');
[wS wC wU wE] = PCA_FIT_FULL_T(WIN,7);
WINY = im2col(Y,[nf 1],'sliding');
WINY = WINY((end-1)/2,:);
SIG = [];
for e = 1:size(wC,1)
    tmp = col2im(wC(e,:),[nf 1],size(X),'sliding');
    %tmp = sum(tmp,1);
    SIG = [SIG;tmp];
end
%% NOPE
close all
netF = patternnet(8);
netF = feedforwardnet([3 8]);
netF.layers{end}.transferFcn = 'logsig';
netF = train(netF,wC,WINY);
fr = sim(netF,wC);
fr = col2im(fr,[nf 1],size(X),'sliding');
frt = col2im(WINY,[nf 1],size(X),'sliding');
%fr = bsxfun(@minus,fr,mean(fr,1));
%% NOPE
close all
tr = [];
for e = 1:size(fr,2)
    os = X(16:(end-15),e);
    [ridx] = find(frt(:,e));
    if ~isempty(ridx)
        tmp = zeros(size(fr,1),1);
        sidx = find(fr(:,e) > graythresh(fr(:,e)));
        sidx = find(fr(:,e) > .6);
        %[~,sidx] = max(fr(:,e));
        tmp(sidx(1)) = 1;
        
        plot(tmp,'r')
        hold on
        plot(os*20000,'b');
        plot(frt(:,e),'k')
        plot(fr(:,e),'g')
        title(KJ{e})
        drawnow
        
        waitforbuttonpress
        
        tr = [tr;[ridx(1) sidx(1)]];
        hold off
    end
end
close all
plot(tr(:,1),tr(:,2),'.');
%% STACKY0
toM = {};
%oM{1} = X2;

toM{1} = bsxfun(@minus,X,mean(X,1));
toM{2} = bsxfun(@minus,X2,mean(X2,1));
toM{3} = bsxfun(@minus,X3,mean(X3,1));
toM{4} = cumsum(toM{1},1);
toM{5} = cumsum(toM{2},1);
toM{6} = cumsum(toM{3},1);
%[~,toM{4}] = gradient(X);
%[~,toM{5}]= gradient(X2);
%[~,toM{6}]= gradient(X3);

%% STACKY2_v2
nf = 11;
nc = 3;
clear wS wC wU wE;
gM = [];
%nc = [5 5 5 5 5 5 4];
nc = [3 3 3 3 3 3 5];
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


%net = feedforwardnet([5 10 8]);
%net = feedforwardnet([20]);
net = patternnet(10);
net.layers{end}.transferFcn = 'logsig';
net = train(net,gM,WINY,'useParallel','yes');
wC{end+1} = sim(net,gM);


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
%% SCORE
[germ frame] = score(toMG,toM,wUG,wEG,wUF,wEF,netG{25},net,nf,nc,ncg);
frame(:,germ==0) = 0;
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
%% TRACK
close all
%[~,TR1] = max(frame,[],1);
%[~,TR2] = max(handD,[],1);
deltaA = abs(tr(:,1) - tr(:,2));
[J,midx] = sort(deltaA,'descend');
viewData(mNM{midx(1)},vecGG(:,:,midx(1))');

%% align genotype SCORE2
for e = 1:numel(subJK)
    geno{e} = GT.get(subJK{e});
    
end
clear para
close all
s
UQ = unique(geno);
for u = 1:numel(UQ)
    %{
    fidx = strcmp(geno,UQ{u});
    subD = algoD(:,fidx);
    binA = mean(subD,2);
    subD = handD(:,fidx);
    binH = mean(subD,2);
    binH = cumsum(binH,1);
    %}
    
    
    
    fidx = strcmp(geno,UQ{u});
    subD = frame(:,fidx);
    ALGO(:,u) = mean(subD,2);
    subD = handD(:,fidx);
    binH = mean(subD,2);
    HAND(:,u) = cumsum(binH,1);
    
    
    sig = ALGO(:,u)';
    xlab = 1:size(ALGO,1);
    [J xval] = min(abs(sig - mean(sig)));
    [para{u}] = fminsearch(@(X)mySigmoid_ver0((1:numel(sig)),X,sig),[sig(end) .5 xlab(xval)]); 
    [e,ALGOf(:,u)] = mySigmoid_ver0(xlab',para{u});
end
CL = {'r' 'g' 'b' 'c' 'k' 'y' 'r--' 'g--' 'b--'};
figure;
hold on
for u = 1:size(HAND,2)
    plot(HAND(:,u),CL{u});
    
end
legend(UQ)
figure
hold on
for u = 1:size(ALGO,2)
    plot(ALGOf(:,u),CL{u});
    
end
legend(UQ)
r = gradient(ALGOf');
[mr,fr] = max(r,[],2);
H = bsxfun(@times,r,sum(r,2).^-1);
E = -sum(log(H).*H,2);
csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/bioM.csv',[ALGOf(end,:)' mr fr E]);
csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/bioD.csv',ALGOf);
csvwrite('/mnt/snapper/nate/mirror_images/maizeData/jgustin/biogD.csv',r);
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
%% complex
finalGerm = hmm_node('finalGerm');
nodeSet = {};

WINY = im2col(Y,[nf 1],'sliding');
WINY = WINY((end-1)/2,:);
WINY = col2im(WINY,[nf 1],size(X),'sliding');
disSet2 = {};
disSet = {};
nodeSet = {};


hmm = my_hmm();

%hmm.addNode(nodeSet{1,1});
%hmm.addNode(nodeSet{2,1});



for s = 1:size(WINY,1)
    UQ = unique(WINY(s,:));
    
    tmpV = squeeze(SIG(s,:,:));
    tmpV2 = squeeze(SIG2(s,:,:));
    
    
    UQ = 1:2;
    for u = 1:numel(UQ)
        nodeSet{u,s} = hmm_node(['stage:' num2str(s) ':state' num2str(UQ(u))]);
        hmm.addNode(nodeSet{u,s});
        if s ~= 1
            fidxT0 = find(all(bsxfun(@eq,WINY((s-1),:),[0]),1));
            fidx0 = find(all(bsxfun(@eq,WINY((s-1):s,:),[0;(UQ(u)-1)]),1));
            zeroToUQ = constantTransitionFunction(numel(fidx0)/numel(fidxT0));

            fidxT1 = find(all(bsxfun(@eq,WINY((s-1),:),[1]),1));
            fidx1 = find(all(bsxfun(@eq,WINY((s-1):s,:),[1;(UQ(u)-1)]),1));
            oneToUQ = constantTransitionFunction(numel(fidx1)/numel(fidxT1));
        
        
            nodeSet{1,s-1}.attachNode(nodeSet{u,s},zeroToUQ);
            nodeSet{2,s-1}.attachNode(nodeSet{u,s},oneToUQ);
        end

        
        preCov = cov(tmpV(WINY(s,:)==(UQ(u)-1),:));
        preU = mean(tmpV(WINY(s,:)==(UQ(u)-1),:),1);
        if any(isnan(preU(:))) | (sum(WINY(s,:)==(UQ(u)-1)) <= 2)
            %preU = mean(tmpV,1);
            %preCov = cov(tmpV);
            preU = zeros(1,size(tmpV,2));
            preCov = eps*eye(2);
        end
        disSet{u,s} = myProb(preU,preCov);
       
        
        preCov2 = cov(tmpV2(WINY(s,:)==(UQ(u)-1),:));
        preU2 = mean(tmpV2(WINY(s,:)==(UQ(u)-1),:));
        if any(isnan(preU2(:))) | (sum(WINY(s,:)==(UQ(u)-1)) <= 2)
            preU2 = zeros(1,size(tmpV,2));
            preCov2 = eps*eye(2);
        end
        disSet2{u,s} = myProb(preU2,preCov2);
        
        
        nodeSet{u,s}.attachDistribution(disSet{u,s},1);
        nodeSet{u,s}.attachDistribution(disSet2{u,s},2);
    end
    
    
end
hmm.dn = [1 1];
%% basic massive
clear postU preU postCov preCov preDis postDis
preIdx = find(WINY==0);
postIdx = find(WINY==1);


postGerm = hmm_node('postGerm');
preGerm = hmm_node('preGerm');


for m = 1:numel(wC)
    
    preCov{1,m} = cov(wC{m}(:,preIdx)');
    preU{1,m} = mean(wC{m}(:,preIdx)');
    postCov{2,m} = cov(wC{m}(:,postIdx)');
    postU{2,m} = mean(wC{m}(:,postIdx)');
    
    preDis{m} = myProb(preU{1,m},preCov{1,m});
    postDis{m} = myProb(postU{2,m},postCov{2,m});
    %{
    preDis{m}.fitToGMM(wC{m}(:,preIdx)',2,1);
    postDis{m}.fitToGMM(wC{m}(:,postIdx)',2,1);
    %}
    
    
    
    preGerm.attachDistribution(preDis{m},m);

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
preTopost = constantTransitionFunction(.5);
postTopost = constantTransitionFunction(1);


preGerm.attachNode(preGerm,myTi);
preGerm.attachNode(postGerm,myT);
postGerm.attachNode(postGerm,postTopost);

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
close all
tr = [];
subK = {};
IDX = [];
n = figure;
cor = figure;
CV = [];
mea = [];
GT = [];
AT = [];
stdM = [];
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
    fidx2 = find(frt(:,e)==1);
    if ~isempty(fidx) & ~isempty(fidx2)
        tr = [tr;[fidx2(1) fidx(1)]];
        subK{end+1} = KJ{e};
        IDX = [IDX ; e];
        
        
        GT = [GT frt(:,e)];
        AT = [AT gidx'==2];
    else
        %break
    end
    figure(n);
    os = X((nf-1)/2:(end-(nf-1)/2),e);
    os2 = X2((nf-1)/2:(end-(nf-1)/2),e);
    os3 = X3((nf-1)/2:(end-(nf-1)/2),e);
    MAG = max([max(os) max(os2)]);
    plot(MAG*frt(:,e),'b');hold all;plot(MAG*(gidx-1),'k');%plot(os*5000 > .05,'g');
    plot(os,'b');
    plot(os2,'r');
    plot(os3,'g');
    if e < 20
        %waitforbuttonpress
    end
    hold off
    drawnow
    e
    CV = [CV corr(tr(:,1),tr(:,2))];
    delta = tr(:,1) - tr(:,2);
    mea = [mea mean(delta)];
    stdM = [stdM std(delta)];
    figure(cor);
    errorbar(mea,stdM);
    drawnow
    
    %waitforbuttonpress
end
figure;
plot(tr(:,1),tr(:,2),'.')
corr(tr(:,1),tr(:,2))
%%
delta = tr(:,1) - tr(:,2);
mean((delta))
[J,midx] = max(abs(delta));
[J,midx] = sort(abs(delta),'descend')
subK{midx(1:10)}
%%

[eS eC eU eE eL eERR eLAM] = PCA_FIT_FULL_T(X,5);
fFrames = 30;
lFrames = 20;
fMean = mean(X(1:fFrames,:));
lMean = mean(X((end-lFrames):end,:));

%[eS SIG eU eE eL eERR eLAM] = PCA_FIT_FULL_T(sort(X,1),3);
[Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(sort(X,1)',Y2',4);

%Yp = [ones(size(X',1),1) sort(X,1)']*beta;
%Yp = Yp';
newD = bsxfun(@minus,sort(X,1)',mean(sort(X,1))');
newD = newD*Weights;
SIG = sort(X,1);
SIG = X;
lambda = mynLDA(SIG',Y2,1,3);

V = bsxfun(@minus,SIG',mean(SIG',1))*lambda;
Yt = zeros(2,numel(Y2));
Yt(2,find(Y2==1)) = 1;
Yt(1,find(Y2==0)) = 1;
%Yt = logical(Yt);
%newD = [];
newD = [newD V];
newD = V;
%newD = [newD abs(lMean - fMean)' V];
%newD = [abs(lMean - fMean)'];
%net = patternnet(20);
clear net cm c
for r = 1:1000
net{r} = feedforwardnet([3 10 10 3]);
net{r}.layers{end}.transferFcn = 'logsig';
net{r} = patternnet(20);
%net{r}.trainFcn = 'trainscg';
%net{r}.performFcn = 'crossentropy';
%eC = SIG;
%newD = [];
%eC = X;
%eC = [];
%newD = [];
%SIG = sort(X,1);
%SIG = max(X,[],1) - min(X,[],1);
%eC = [min(X,[],1);max(X,[],1)];
%eC = eC.*flipud(eC).^-1;
net{r} = train(net{r},[eC;newD'],Yt);
%net.train([eC;newD'],Yt);

Yp = sim(net{r},[eC;newD']);


%Yp = round(Yp);
Yp2 = zeros(si
ze(Yp));
Yp2(2,find(Yp==1)) = 1;
Yp2(1,find(~(logical(Yp)))) = 1;
%
[c(r),cm{r},ind,per] = confusion(Yp,Yt);
cm{r}
min(c)
end
[J,midx] = min(c);
gNet = net{midx};
%%
    
    for f = 1:numel(d)
        if ~isempty(strfind(em{e,iplantName},d(f).experimentName))
            try
                arrayName = em{e,posName};
                arrayName = ['; ' arrayName ';'];
                fidx = strfind(arrayName,';');
                for l = 1:(numel(fidx)-1)
                    snip = arrayName((fidx(l)+1):(fidx(l+1)-1));
                    sidx = strfind(snip,':');
                    letterRange1 = snip(2);
                    letterRange2 = snip(sidx(1)+1);
                    numberRange1 = str2num(snip(3:(sidx(1)-1)));
                    numberRange2 = str2num(snip((sidx(1)+2):end));
                    gidx = strfind(em{e,iplantName},filesep);
                    cn = numberRange1:numberRange2;
                    culN = em{e,6};
                    for rn = 1:numel(LL)
                        for c = 1:numel(cn)
                            %key = [em{e,iplantName}((gidx(end)+1):end) '-' lower(em{e,hcName}(1)) '-' LL{rn} num2str(cn(c))];
                            %flag  = ~strcmp(ml{e,6},'15A-6197-02');

                            key = [em{e,iplantName}((gidx(end)+1):end) '-' lower(em{e,hcName}(1)) '-' LL{rn} num2str(cn(c))];
                            
                            
                            
                            key = lower(key);
                            key;

                            
                            value = LT.get(key);
                            em{e,FRAME_SCORE};
                            if ~isnan(em{e,FRAME_SCORE})
                                tmp = zeros(size(value));
                                %waitforbuttonpress
                                tmp(em{e,FRAME_SCORE}) = 1*max(value);
                                plot(value);
                                hold on
                                plot(tmp);
                                hold off
                                X = [X value];
                                Y = [Y tmp];
                                drawnow
                                %waitforbuttonpress
                            end
                            %{
                            deltaT = find(value);
                            genoTmp = ml{e,genoName};
                            %posTmp = EM.get(key);
                            vec = GT.get(genoTmp);
                            vec = [vec value(1:250)];
                            if flag
                                GT.put(genoTmp,vec);
                            end


                            if isempty(deltaT)
                                deltaT = inf;
                            end
                            %corrV = [posTmp,deltaT(1)];
                            %if flag
                            %    MT.put(key,corrV);
                            %end

                            vec = GV.get(genoTmp);
                            %vec = [vec corrV'];
                            %if flag
                            %    GV.put(genoTmp,vec);
                            %end

                            vec = SM2.get(culN);
                            vec = [vec;deltaT(1)];
                            if flag
                                SM2.put(culN,vec);
                            end
                            %}
                        end
                    end
                end
            catch ME
                ME
            end
            %lower(ml{e,hcName}(1));
        end
    end
end
%% read data for emergnece
clear d
d1 = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/20170131_Camera3.csv');
d2 = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/20170131_Camera4.csv');

d(1).experimentName = '20170131_Camera3';
d(2).experimentName = '20170131_Camera4';
d(1).data = cell2mat(d1);
d(2).data = cell2mat(d2);
%% emergnece master list
ml = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/EmergenceAssay_Master_list.csv');
%% auto generate d structe
FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
for e = 1:numel(FileList)
    fileName = FileList{e};
    fidx = strfind(fileName,'_');
    [pth,fn,ext] = fileparts(FileList{e});
    
    
    tmp = readtext(FileList{e});
    
    d(e).data = cell2mat(tmp);
    d(e).experimentName = fn;
end
%% other data

em = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/NAM_parents_ear_map_xyCoords(1).csv');
%% for swelling

FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/kernel_swellData/return/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
kp = [];
pkp = [];
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'swell--'))
        kp = [kp e];
    end
     if ~isempty(strfind(FileList{e},'para--'))
        pkp = [pkp e];
    end
end
pFileList = FileList(pkp);
FileList = FileList(kp);

close all
im = readtext('/home/nate/Downloads/imbibition_assay_master_list(1).csv');
LLI = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K'};
import java.util.Map.*;
import java.util.HashMap;
SM = HashMap();
for e = 2:size(im,1)
    earKey = im{e,6}
    scannerKey = im{e,4};
    arrayName = im{e,10};
    arrayName = ['; ' arrayName ';'];
    fidx = strfind(arrayName,';');
    fidx = strfind(arrayName,';');
    for l = 1:(numel(fidx)-1)
        snip = arrayName((fidx(l)+1):(fidx(l+1)-1));
        sidx = strfind(snip,':');
        letterRange1 = snip(2);
        letterRange2 = snip(sidx(1)+1);
        numberRange1 = str2num(snip(3:(sidx(1)-1)));
        numberRange2 = str2num(snip((sidx(1)+2):end));
        
        rowN = find(strcmp(LLI,letterRange1));
        for l = 1:numel(FileList)
            if ~isempty(strfind(FileList{l},['--' num2str(rowN)])) & ~isempty(strfind(FileList{l},scannerKey))
                sd = csvread(FileList{l});
                pd = csvread(pFileList{l});
                [keep] = manfredFilter(sd');
                plot(sd,'k')
                hold on
                plot(sd(:,keep),'r')
                hold off
                drawnow
                %waitforbuttonpress
                earMean = mean(sd(:,keep),2);
                earMean = mean(pd(keep,2));
                SM.put(earKey,earMean);
            end
        end
    end
end
%% parse ear map
cameraName = 8;
posML2 = 9;
import java.util.Map.*;
import java.util.HashMap;
EM = HashMap();
for e = 2:size(em,1)
    if ~strcmp(em{e,cameraName},'NA')
        key = ['20170131_' lower(em{e,cameraName}) '-' lower(em{e,posML2}(1)) '-' lower(em{e,posML2}(2:end))];
        value = em{e,12} - em{e,14};
        LEN = em{e,15} - em{e,14};
        value = value/LEN;
        if strcmp(key,'20170131_camera4-d-g12')
            value
        end
        if isempty(value)
            key;
        end
        key;
        EM.put(key,value);
    end
end
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
LT = HashMap();
FNK = HashMap();
for e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    fidx = strfind(FileList{e},'_');
    gidx = strfind(FileList{e},'.');
    tmpK = FileList{e}((fidx(end)+1):(gidx(end)-1));
    newN = trK.get(tmpK);
    newDataKey{e} = strrep(nm,tmpK,newN);
    fidx = strfind(newDataKey{e},'_');
    newDataKey{e}(fidx(end)) = '-';
    
    
    tmpD = O{e}';
    for k = 1:size(tmpD,2)
        tmpD(:,k) = imfilter(tmpD(:,k),fspecial('average',[11 1]),'replicate');
        %tmpD(:,k) = subBL(tmpD(:,k),50);
    end
    
    LT.put(lower(newDataKey{e}),tmpD);
    FNK.put(lower(newDataKey{e}),FileList{e});
end
%% pull data from old measures
K = {};
for e = 1:numel(d)
    for f = 1:size(d(e).ampSig,2)
        key = lower([d(e).experimentName '-' LABELS{f}]);
        [dL var1] = subBL(d(e).ampSig(:,f),50);
        dL = imfilter(dL,fspecial('average',[11 1]),'replicate');
        dL = dL / sum(dL);
        %dL = dL / max(dL);
        [ampSig  var2]= subBL(d(e).redSig(:,f),50);
        ampSig = imfilter(ampSig,fspecial('average',[11 1]),'replicate');
        ampSig = ampSig / sum(ampSig);
        %ampSig = ampSig / max(ampSig);
        %other = .5*(dL+ampSig);
        other = dL.*ampSig;
        K{end+1} = key;
        value = [dL ampSig other];
        %value = [dL;ampSig];
        LT.put(key,value);
        
    end
end
done = 1
%% align data from mastersheet
iplantName = 27+1;
hcName = 20+1;
posName = 10+1;
genoName = 7+1;
import java.util.Map.*;
import java.util.HashMap;
GT = HashMap();
MT = HashMap();
GV = HashMap();
SM2 = HashMap();
for e = 2:size(ml,1)
    for f = 1:numel(d)
        if ~isempty(strfind(ml{e,iplantName},d(f).experimentName))
            
            arrayName = ml{e,posName};
            arrayName = ['; ' arrayName ';'];
            fidx = strfind(arrayName,';');
            for l = 1:(numel(fidx)-1)
                snip = arrayName((fidx(l)+1):(fidx(l+1)-1));
                sidx = strfind(snip,':');
                letterRange1 = snip(2);
                letterRange2 = snip(sidx(1)+1);
                numberRange1 = str2num(snip(3:(sidx(1)-1)));
                numberRange2 = str2num(snip((sidx(1)+2):end));
                gidx = strfind(ml{e,iplantName},filesep);
                cn = numberRange1:numberRange2;
                culN = ml{e,6};
                for rn = 1:numel(LL)
                    for c = 1:numel(cn)
                        key = [ml{e,iplantName}((gidx(end)+1):end) '-' lower(ml{e,hcName}(1)) '-' LL{rn} num2str(cn(c))];
                        flag  = ~strcmp(ml{e,6},'15A-6197-02');
                        
                        key = lower(key);
                        key
                        value = LT.get(key);
                        deltaT = find(value);
                        genoTmp = ml{e,genoName};
                        %posTmp = EM.get(key);
                        vec = GT.get(genoTmp);
                        vec = [vec value(1:250)];
                        if flag
                            GT.put(genoTmp,vec);
                        end
                        
                        
                        if isempty(deltaT)
                            deltaT = inf;
                        end
                        %corrV = [posTmp,deltaT(1)];
                        %if flag
                        %    MT.put(key,corrV);
                        %end
                        
                        vec = GV.get(genoTmp);
                        %vec = [vec corrV'];
                        %if flag
                        %    GV.put(genoTmp,vec);
                        %end
                        
                        vec = SM2.get(culN);
                        vec = [vec;deltaT(1)];
                        if flag
                            SM2.put(culN,vec);
                        end
                        
                    end
                end
            end
            lower(ml{e,hcName}(1));
        end
    end
end
%%
close all
keys = MT.keySet;
itr = keys.iterator();
STACK = [];

while itr.hasNext()
    key = itr.next();
    value = MT.get(key);
    
    if ~isinf(value(2))
        plot(value(1),value(2),'.')
        STACK = [STACK;value'];
        hold on
    end
end
%%
close all
keys = SM.keySet;
itr = keys.iterator();
while itr.hasNext()
    key = itr.next();
    val1 = SM.get(key);
    val2 = SM2.get(key);
    val2(isinf(val2)) = [];
    plot(mean(val2),val1,'.')
    hold on
end
%%
close all
keys = GT.keySet;
itr = keys.iterator();
%h1 = figure;
%h2 = figure;
cnt = 1;
mPDF = [];
iPDF = [];
while itr.hasNext()
    %for e = 1
    %key = itr.next();
    key = itr.next();
    vec = GT.get(key);
    %{
    %whole = GV.get(key);
    
    %nonGerm = isinf(whole(2,:));
    %figure(h2);
    %sum(nonGerm)
    f = ksdensity(whole(1,nonGerm),linspace(0,1,1000));
    f = f / sum(f);
    plot(linspace(0,1,1000),f)
    hold all
    plot(whole(1,:),ones(1,numel(nonGerm))*cnt*.001,'.');
    %}
    
    
    
    
    
    %{
    figure(h1);
    
    
    
    
    whole(:,isinf(whole(2,:))) = [];
    plot(whole(1,:),whole(2,:),'.');
    hold all
    %}
    
    no_germN = sum(all(vec==0,1));
    germN = size(vec,2) - no_germN;
    germP(cnt,:) = [germN no_germN]/size(vec,2);
    
    ipdf = gradient(imfilter(mean(vec,2),fspecial('average',[5 1]),'replicate'));
    pdf = imfilter(mean(vec,2),fspecial('average',[5 1]),'replicate')
    ipdf = ipdf / 30  * 60 ;
    ipdf(ipdf < 0) = 0;
    xlab = 1:numel(ipdf);
    xlab = (xlab*30/60 + 80)/24;
    
    
    [J xval] = min(abs(pdf - mean(pdf)));
    [para{cnt}] = fminsearch(@(X)mySigmoid_ver0(xlab',X,pdf),[pdf(end) .5 xlab(xval)]); 
    [e,yp] = mySigmoid_ver0(xlab',para{cnt});
    ypg = gradient(yp);
    %{
    type = 'logit';
    type = 'probit';
    n = size(vec,2)*ones(size(pdf));
    xlab = 1:size(vec,1);
    toFit = sum(vec,2);
    fidx = find(toFit);
    [logitCoef,dev] = glmfit(xlab(fidx(1):end)',[toFit(fidx(1):end) n(fidx(1):end)],'binomial','link',type);
    logitFit = glmval(logitCoef,xlab',type,'size', n);
    plot(xlab,logitFit);
    %}
    
    
    figure;
    [AX,H1,H2] = plotyy(xlab,ipdf*100,xlab,pdf);
    hold(AX(2))
    hold(AX(1))
    plot(AX(2),xlab',yp,'r')
    plot(AX(1),xlab',ypg*100/ 30  * 60 ,'r')
    set(AX(2),'YTick',linspace(0,1,11));
    axis(AX(2),[3 9 0 1]);
    axis(AX(1),[3 9 0 10]);
    set(AX(1),'YTick',linspace(0,10,11));
    set(get(AX(1),'Ylabel'),'String','percent germ per hour') 
    set(get(AX(2),'Ylabel'),'String','total percent germ')
    xlabel('days after planting');
    
    LEG{cnt} = key;
    mPDF(:,cnt) = yp;
    iPDF(:,cnt) = ypg;
    %title([key '-' num2str(size(vec,2))])
    title([key])
    %{
     hold on
    figure;
    plot(xlab',yp)
    hold on
    plot(xlab',pdf,'r')    
    %}
    cnt = cnt + 1;
end
%%
CL = {'r' 'g' 'b' 'm' 'c' 'k' 'r--' 'b--' 'g--'};
h = figure;
ax(1) = axes();
ax(2) = axes('YAxisLocation','right');
for e = 1:size(mPDF,2)
    plot(ax(1),xlab',mPDF(:,e),CL{e});
    hold(ax(1),'on');
    plot(ax(2),xlab',100*iPDF(:,e)/ 30  * 60 ,CL{e});
    hold(ax(2),'on');
end

axis(ax(1),[3 9 0 1]);
set(ax(1),'YTick',linspace(0,1,11));
axis(ax(2),[3 9 0 10]);
set(ax(2),'YTick',linspace(0,10,11));
set(ax(2),'YAxisLocation','right');
set(ax(2),'Color','none');
ylabel(ax(1),'Percent Germination')
ylabel(ax(2),'Percent/hour Germination')
xlabel(ax(1),'Time (days)')
legend(LEG)
%%

figure;
bar(germP,'stack')



