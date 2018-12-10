preData = readtext('/home/nate/Downloads/ScoreMergeEML.csv');
%%
preData = readtable('/home/nate/Downloads/ScoreMergeEML(1).csv');
%%
unique(preData.experiment_code)
%% split out the project via experiemtn code
pidx = find(strcmp(preData.experiment_code,'Mo17 Cold Senstivity Test'));
%%
pidx = find(strcmp(preData.experiment_code,'Wisconsin_Standardization_Test'));

%%
ecode = preData.experiment_code(pidx);
scores = preData.EmergenceHour_MachineScore_(pidx);
genotype = preData.Genotype(pidx);
isolate = preData.Isolate(pidx);
cultureNumber = preData.Culture_Number(pidx);
treatment =preData.Seed_treatment_1(pidx);
season = {};
for e = 1:numel(cultureNumber)
    season{e} = cultureNumber{e}(1:3);
end
TIME = linspace(0,max(scores(~isinf(scores))),1000);
step = zeros(numel(scores),numel(TIME));
for e = 1:numel(scores)
    if ~isinf(scores(e))
        [~,idx] = min(abs(TIME-scores(e)));

        shape = normcdf(TIME,TIME(idx),4);
        step(e,idx:end) = 1;
        step(e,:) = shape;
        
    end
end
%%
close all
ukey = {};
for e = 1:numel(season)
    ukey{e} = [season{e} '_' isolate{e} '_' treatment{e}];
end
UQ = unique(ukey);
rm = [];
for u = 1:numel(UQ)
    if ~isempty(strfind(UQ{u},'5C'))
        rm(u) = true;
    end
end
UQ(find(rm)) = [];
uS = [];
for u = 1:numel(UQ)
    uidx = strcmp(ukey,UQ{u});
    uS(u,:) = mean(step(uidx,:),1);
end
plot(TIME/24,uS');
%%
groups = 1:size(uS,1);
%%
close all
plot(TIME,uS(3,:))
%%

D = uS(1:2,:);
D = uS(1:5,:);
D = uS;
% time for curves
TM = repmat(TIME,[size(D,1) 1]);

groups = (1:size(D,1))';
groups = repmat(groups,[1 size(D,2)]);
options = statset('Display','iter');
ridx = isnan(D(:));
TM(ridx) = [];
groups(ridx) = [];
D(ridx) = [];
finalGuess = mean(D(:,end));

funcToFit = @(X,tm)X(1)*(1+exp(-X(2)*(tm-X(3)))).^-1;
MX = max(mean(D,1));
[~,xo] = min(abs(mean(D,1) - MX/2));

%[fixed,PSI1,stats,random] = nlmefit(TM(:),D(:),groups(:),[],funcToFit,[MX .1 xo],'Options',options);  
clear delta
for e = 1:size(D,1)
    [beta{e},R{e},J{e},CovB{e},MSE{e},ErrorModelInfo{e}] = nlinfit(TIME,D(e,:),funcToFit,[MX .1 xo]);
    [Ypred(e,:),delta(e,:)] = nlpredci(funcToFit,TIME,beta{e},R{e},'Covar',CovB{e});
end
%%

%% single
for e = 1:numel(UQ)
    fidx = find(strcmp(UQ{e},ukey));
    D = uS(1:2,:);
    D = uS(1:5,:);
    D = uS;
    D = step(fidx,:);
    ridx = all(D == 0,2);
    MX(e) =  mean(~ridx);
    D(ridx,:) = [];
    sTIME = TIME;
    oD = D;
    D = D(:,1:10:end);
    sTIME = sTIME(:,1:10:end);
    % time for curves
    TM = repmat(sTIME,[size(D,1) 1]);

    groups = (1:size(D,1))';
    groups = repmat(groups,[1 size(D,2)]);
    options = statset('Display','iter');
    ridx = isnan(D(:));
    TM(ridx) = [];
    groups(ridx) = [];
    D(ridx) = [];
    finalGuess = mean(D(:,end));

    funcToFit = @(X,tm)X(1)*(1+exp(-X(2)*(tm-X(3)))).^-1;
    funcToFit = @(X,tm)normcdf(tm,X(1),4);
    %MX = max(mean(D,1));
    [~,xo] = min(abs(mean(D,1) - MX(e)/2));

    [fixed{e},PSI1{e},stats{e},random{e}] = nlmefit(TM(:),D(:),groups(:),[],funcToFit,[xo],'Options',options);
 

    %{


    close all

    Ypred = MX(e)*funcToFit(fixed{e},TIME);
    qq = MX(e)*funcToFit(std(random{e},1,2)+fixed{e},TIME);
    qn = MX(e)*funcToFit(-std(random{e},1,2)+fixed{e},TIME);

    

    plot(Ypred,'r')
    hold on
    plot(qq,'r--');
    plot(qn,'r--');
    hold on
    plot(uS(e,:),'k')
    hold off
   % waitforbuttonpress
    drawnow
    %}

    for m = 1:size(random{e},2)
        plot(oD(m,:),'k')
        hold on
        yi = funcToFit(fixed{e}+random{e}(m),TIME);
        plot(yi,'r')

    end


end
%%
for e = 1:numel(fixed)
    close all

    Ypred = MX(e)*funcToFit(fixed{e},TIME);
    qq = MX(e)*funcToFit(std(random{e},1,2)+fixed{e},TIME);
    qn = MX(e)*funcToFit(-std(random{e},1,2)+fixed{e},TIME);

    

    plot(Ypred,'r')
    hold on
    plot(qq,'r--');
    plot(qn,'r--');
    hold on
    plot(uS(e,:),'k')
    hold off
    waitforbuttonpress
end
%% 
[q1 q2] = nlpredci(funcToFit,mean(TM),fixed,PSI1,'Jacobian',stats);
funcToFit(fixed,mean(TM))
%%
emergenceTable = table();
TH = .8;
dT = mean(gradient(TIME));
for e = 1:size(fitF,1)
    df = gradient(fitF(e,:));
    [mx wmx] = max(df);
    emergenceTable{e,'percentEmergence'} = fitF(e,end);
    emergenceTable{e,'maxPercent'} = mx;
    emergenceTable{e,'whenMax'} = wmx;
    sdf = sort(df,'descend');
    [~,midx] = min(abs(cumsum(sdf) - TH*fitF(e,end)));
    segment = df > sdf(midx);
    emergenceTable{e,'duration'} = sum(segment)*dT;
    [AX,H1,H2] = plotyy(TIME,fitF(e,:),TIME,df);
    figure;
    plot(AX,uS(e,:),'k');
end
plot(df)
hold on
plot(segment)
%%
emergenceTable = table();
TH = .8;
dT = mean(gradient(TIME));
for e = 1:numel(UQ)
    uidx = find(strcmp(UQ{e},ukey));
    subD = scores(uidx);
    subC = cultureNumber(uidx);
    UQ2 = unique(subC);
    tempy = [];
    for u2 = 1:numel(UQ2)
        sidx1 = strcmp(UQ2{u2},subC);
        subby = subD(sidx1);
        tempy(u2) = mean(~isinf(subby));
    end
    %df = gradient(fitF(e,:));
    %[mx wmx] = max(df);
    emergenceTable{UQ{e},'percentEmergence'} = mean(tempy);
    emergenceTable{UQ{e},'percentEmergenceERR'} = std(tempy);
    emergenceTable{UQ{e},'trialSize'} = size(subD,1);
    %emergenceTable{e,'maxPercent'} = mx;
    emergenceTable{UQ{e},'whenMax'} = mean(subD(~isinf(subD)));
    %{
    sdf = sort(df,'descend');
    [~,midx] = min(abs(cumsum(sdf) - TH*fitF(e,end)));
    segment = df > sdf(midx);
    %}
    emergenceTable{UQ{e},'duration'} = std(subD(~isinf(subD)));



    %[AX,H1,H2] = plotyy(TIME,fitF(e,:),TIME,df);
    %figure;
    %plot(AX,uS(e,:),'k');
end

plot(df)
hold on
plot(segment)
%%
close all
P = funcToFit([.9 .1 400],TIME);
plot(P)
%%




