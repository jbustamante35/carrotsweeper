clear all
%%
% read in file from bessie
%L = readtext('/mnt/snapper/Bessie/NMDAreceptor/NMDA CSV FILES/06feb18length.csv');
L = readtext('/mnt/snapper/Bessie/NMDAreceptor/NMDA CSV FILES/Length.Compiled.csv');
T = readtable('/mnt/snapper/Bessie/NMDAreceptor/NMDA CSV FILES/Length.Compiled.csv');
%%
L = readtext('/mnt/snapper/Bessie/NMDAreceptor/NMDA CSV FILES/26jan18tipangle.csv');

%% extract treatment and genotype and p as phenotype
treatment = L(:,2);
genotype = L(:,3);
p = L(:,4:end);
%% make labels for treatment computer consistant
types = {'glu','gly','min'};
clear TYPE
for e = 1:numel(treatment)
    for t = 1:numel(types)
        TYPE(e,t) = ~isempty(strfind(treatment{e},types{t}));
    end
end


% replace glu.gly
for e = 1:numel(treatment)
    if strcmp(treatment{e},'glu.gly')
        treatment{e} = '1mmglu.gly';
    end
    if strcmp(treatment{e},'min')
        treatment{e} = '0mmmin';
    end
end

% amount
for e = 1:numel(treatment)
    if ~isempty(strfind(treatment{e},'0mm'))
        amount{e} = '0mm';
    elseif ~isempty(strfind(treatment{e},'1mm'))
        amount{e} = '1mm';
    elseif ~isempty(strfind(treatment{e},'100um'))
        amount{e} = '100um';
    else
        'no'
    end
end


% phenotype
for e = 1:numel(p)
    if isempty(p{e})
        p{e} = NaN;
    end
end
phenotype = cell2mat(p);

%% clean out missing data
mask = find(any(phenotype == 0 | isnan(phenotype),2) | any(abs(diff(phenotype,1,2)) > 5,2));
%mask = find(any(phenotype == 0 | isnan(phenotype),2) | any(abs(diff(phenotype,1,2)) > .5,2));
genotype(mask) = [];
amount(mask) = [];
treatment(mask) = [];
TYPE(mask,:) = [];
phenotype(mask,:) = [];
mask = phenotype == 0 | isnan(phenotype);
any(mask,1);


size(genotype)
size(treatment)
size(phenotype)

phenotype = bsxfun(@minus,phenotype,phenotype(:,1));
%% PCA on phenotype
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(phenotype,3);
%% obtain the average rate
rate = mean(diff(phenotype,1,2),2);
SIM = rate*(1:size(phenotype,2));
DELTA = phenotype - SIM;
sDELTA = imfilter(DELTA,fspecial('average',[1 21]),'replicate');
[dS dC dU dE dL dERR dLAM] = PCA_FIT_FULL(sDELTA,1);
%% parse treatment
phenotypeData = table;
for e = 1:numel(treatment)
    
    
    if strfind(genotype{e},'wt')
        phenotypeData{e,'genotype'} = {'wt'};
    else
        phenotypeData{e,'genotype'} = {'NMDA'};
    end
    
    if strfind(genotype{e},'2')
        phenotypeData{e,'grow_out'} = {'2'};
    else
        phenotypeData{e,'grow_out'} = {'1'};
    end
    
    if strfind(treatment{e},'glu')
        if ~isempty(strfind(amount{e},'mm'))
            phenotypeData{e,'Glu_Conc_um'} = 1000;
        elseif ~isempty(strfind(amount{e},'um'))
            phenotypeData{e,'Glu_Conc_um'} = 100;
        else
            phenotypeData{e,'Glu_Conc_um'} = 0;
        end
    end
    
    if strfind(treatment{e},'gly')
        if ~isempty(strfind(amount{e},'mm'))
            phenotypeData{e,'Gly_Conc_um'} = 1000;
        elseif ~isempty(strfind(amount{e},'um'))
            phenotypeData{e,'Gly_Conc_um'} = 100;
        else
            phenotypeData{e,'Gly_Conc_um'} = 0;
        end
    end
    
    for pc = 1:size(pC,2)
        phenotypeData{e,['PC_' num2str(pc)]} = pC(e,pc);
    end
    
    
    phenotypeData{e,'averageRate'} = rate(e);
    
    phenotypeData{e,'delta'} = dC(e);
    
end
%%
IDX = find(phenotypeData.Glu_Conc_um == 100 & phenotypeData.Gly_Conc_um == 0);
%IDX = find(phenotypeData.Gly_Conc_um == 1000);
%IDX = find(phenotypeData.Glu_Conc_um == 1000);
%IDX = find(phenotypeData.Gly_Conc_um == 100);
%IDX = find(phenotypeData.Glu_Conc_um == 100);
K = phenotypeData.grow_out(IDX);
kp = [];
for e = 1:numel(K)
    kp(e) = strcmp(K{e},'1');
end
kp = find(kp);
GR = []
P = phenotypeData.averageRate(IDX(kp));
G = phenotypeData.genotype(IDX(kp));
for e = 1:numel(G)
    if strcmp(G{e},'wt')
        GR(e) = 1;
    end
end
close all
ksdensity(P(GR==1))
waitforbuttonpress
hold on
ksdensity(P(GR==0))
[a b] = ttest2(P(GR==1),P(GR==0))
%%

model{1} = fitlme(phenotypeData,...
'averageRate ~ 1  + Glu_Conc_um  + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

model{2} = fitlme(phenotypeData,...
'averageRate ~ 1  + Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

model{3} = fitlme(phenotypeData,...
'averageRate ~ 1  + Gly_Conc_um + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

model{4} = fitlme(phenotypeData,...
'averageRate ~ 1  + Glu_Conc_um + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

model{5} = fitlme(phenotypeData,...
'averageRate ~ 1  + Gly_Conc_um + Glu_Conc_um + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

model{6} = fitlme(phenotypeData,...
'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

model{7} = fitlme(phenotypeData,...
'averageRate ~ 1  + (Gly_Conc_um + Glu_Conc_um + Gly_Conc_um:Glu_Conc_um |grow_out:genotype)',...
'FitMethod','ML');


model{8} = fitlme(phenotypeData,...
'averageRate ~ 1  +  Gly_Conc_um + Glu_Conc_um:Gly_Conc_um + (genotype|grow_out)',...
'FitMethod','ML');

model{9} = fitlme(phenotypeData,...
'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (genotype|grow_out)',...
'FitMethod','ML');

model{10} = fitlme(phenotypeData,...
'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (Glu_Conc_um:Gly_Conc_um|genotype) + (genotype|grow_out)',...
'FitMethod','ML');


model{11} = fitlme(phenotypeData,...
'averageRate ~ 1  + (Glu_Conc_um:Gly_Conc_um|genotype:grow_out)',...
'FitMethod','ML');


model{12} = fitlme(phenotypeData,...
'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|genotype:grow_out)',...
'FitMethod','ML');
%% sig models -fixed effects
model{1} = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (genotype|grow_out)','FitMethod','ML');
model{2} = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um + (genotype|grow_out)','FitMethod','ML');
model{3} = fitlme(phenotypeData,'averageRate ~ 1  + Gly_Conc_um + (genotype|grow_out)','FitMethod','ML');
model{4} = fitlme(phenotypeData,'averageRate ~ 1  + Gly_Conc_um + Glu_Conc_um:Gly_Conc_um + (genotype|grow_out)','FitMethod','ML');
options = statset('UseParallel',true);
modelCompare = compare(model{3},model{1},'nsim',1000,'Options',options);
%% pref B
A = fitlme(phenotypeData,'averageRate ~ 1  + (Glu_Conc_um:Gly_Conc_um|grow_out:genotype)','FitMethod','ML');
B = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|grow_out:genotype)','FitMethod','ML');
[~,~,stats] = randomEffects(B);
modelCompare = compare(A,B,'nsim',1000,'Options',options);
%% pref B
A = fitlme(phenotypeData,'averageRate ~ 1  + (Glu_Conc_um:Gly_Conc_um:genotype|grow_out)','FitMethod','ML');
B = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (genotype|grow_out)','FitMethod','ML');
[~,~,stats] = randomEffects(B);
modelCompare = compare(A,B,'nsim',100,'Options',options)
%% perf Ap
Ap = fitlme(phenotypeData,'averageRate ~ 1  + Gly_Conc_um + (Glu_Conc_um:Gly_Conc_um|grow_out:genotype)','FitMethod','ML');
Bp = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (Gly_Conc_um|grow_out:genotype)','FitMethod','ML');
[~,~,stats] = randomEffects(Bp);
modelCompare = compare(Ap,Bp,'nsim',1000,'Options',options);
%% compare B to Ap = B
modelCompare = compare(B,Ap,'nsim',1000,'Options',options);
%% try = App
App = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|grow_out:genotype)','FitMethod','ML');
Bpp = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (Gly_Conc_um|grow_out:genotype)','FitMethod','ML');
[~,~,stats] = randomEffects(Bpp);
modelCompare = compare(App,Bpp,'nsim',100,'Options',options);
%% try = App
App = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|grow_out:genotype)','FitMethod','ML');
Bpp = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (Glu_Conc_um|grow_out:genotype)','FitMethod','ML');
[~,~,stats] = randomEffects(Bpp);
modelCompare = compare(App,Bpp,'nsim',100,'Options',options);
%% try = App
App = fitlme(phenotypeData,'delta ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|grow_out:genotype)','FitMethod','ML');
Bpp = fitlme(phenotypeData,'delta ~ 1  + Glu_Conc_um:Gly_Conc_um + (Gly_Conc_um|grow_out:genotype)','FitMethod','ML');
[~,~,stats] = randomEffects(Bpp);
modelCompare = compare(App,Bpp,'nsim',100,'Options',options);

%% try = App
App = fitlme(phenotypeData,'delta ~ 1 + Glu_Conc_um:Gly_Conc_um + Gly_Conc_um + Glu_Conc_um + (1|grow_out:genotype)','FitMethod','ML');
Bpp = fitlme(phenotypeData,'delta ~ 1 + Glu_Conc_um:Gly_Conc_um + Gly_Conc_um + (Glu_Conc_um|grow_out:genotype)','FitMethod','ML');
[~,~,stats] = randomEffects(Bpp);
modelCompare = compare(App,Bpp,'nsim',1000,'Options',options);
%%
[B,Bnames,statsE] = randomEffects(model{11})
[psi,mse,statsC] = covarianceParameters(model{10})
%%
BEST = fitlme(phenotypeData,'averageRate ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|grow_out:genotype)','FitMethod','ML');
[B,Bnames,stats] = randomEffects(BEST)
%%
tmp = [];
figure;
CL = {'r' 'g' 'b'};
cL = [0 100 1000];
gQ = 'NMDA';
for e = 1:size(phenotypeData,1)
    g = phenotypeData{e,'genotype'};
    g = g{1};
    c = phenotypeData{e,'Glu_Conc_um'};
   
    if strcmp(g,gQ) & c == 0
        g
        c
        
        
        c = phenotypeData{e,'Gly_Conc_um'};
        
        plot(phenotype(e,:),CL{cL==c})
        tmp = [tmp;[phenotypeData.averageRate(e) c]];
        %phenotypeData(e)
    end
   
    hold on
end
figure;

plot(tmp(:,2),tmp(:,1),'.')
axis([-10 1200 0 .8]);
title(gQ);
%%

modelD{1} = fitlme(phenotypeData,...
'delta ~ 1  + Glu_Conc_um  + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

modelD{2} = fitlme(phenotypeData,...
'delta ~ 1  + Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

modelD{3} = fitlme(phenotypeData,...
'delta ~ 1  + Gly_Conc_um + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

modelD{4} = fitlme(phenotypeData,...
'delta ~ 1  + Glu_Conc_um + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

modelD{5} = fitlme(phenotypeData,...
'delta ~ 1  + Gly_Conc_um + Glu_Conc_um + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

modelD{6} = fitlme(phenotypeData,...
'delta ~ 1  + Glu_Conc_um:Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');

modelD{7} = fitlme(phenotypeData,...
'delta ~ 1  +  Gly_Conc_um + Glu_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');


%%
%{
complexModel3 = fitlme(phenotypeData,...
'PC_1 ~ 1  + Glu_Conc_um + Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'FitMethod','ML');
%}
options = statset('UseParallel',true);
for e = 8%1:numel(model)
    for k = 9%1:numel(model)
        modelCompare{e,k} = compare(model{e},model{k},'nsim',1000,'Options',options);
    end
end
%%
for e = 1:numel(model)
    for k = 1:numel(model)
        tmp = double(modelCompare{e,k}.pValue);
        PV(e,k) = tmp(2);
    end
end

%%
phenotypeDataNewWT = table;
CON = linspace(0,1000,5);
genoType = 'wt';
growout = '1';
Gly_Conc_um = 1000;
[GLU,GLY] = ndgrid(linspace(0,5000,50),linspace(0,5000,50));
CON = [GLU(:),GLY(:)];
for e = 1:size(CON,1)
    phenotypeDataNewWT{e,'genotype'} = {genoType};
    phenotypeDataNewWT{e,'grow_out'} = {growout};
    phenotypeDataNewWT{e,'Gly_Conc_um'} = CON(e,2);
    phenotypeDataNewWT{e,'Glu_Conc_um'} = CON(e,1);
end
%%
phenotypeDataNewNMDA = table;
CON = linspace(0,1000,5);
genoType = 'NMDA';
growout = '1';
Gly_Conc_um = 1000;
[GLU,GLY] = ndgrid(linspace(0,5000,50),linspace(0,5000,50));
CON = [GLU(:),GLY(:)];
for e = 1:size(CON,1)
    phenotypeDataNewNMDA{e,'genotype'} = {genoType};
    phenotypeDataNewNMDA{e,'grow_out'} = {growout};
    phenotypeDataNewNMDA{e,'Gly_Conc_um'} = CON(e,2);
    phenotypeDataNewNMDA{e,'Glu_Conc_um'} = CON(e,1);
end
%%
toUSE = 11;
predSlopeWT = predict(model{toUSE},phenotypeDataNewWT);
predSlopeNMDA = predict(model{toUSE},phenotypeDataNewNMDA);
predSlopeWT = reshape(predSlopeWT,size(GLU));
predSlopeNMDA = reshape(predSlopeNMDA,size(GLU));
close all
mesh(linspace(0,5000,50)/1000,linspace(0,5000,50)/1000,predSlopeWT.*(predSlopeWT > 0));
hold on
mesh(linspace(0,5000,50)/1000,linspace(0,5000,50)/1000,predSlopeNMDA.*(predSlopeNMDA > 0));
xlabel('GLU')
ylabel('GLY')

%%
toUSE = 9;
predSlopeWT = predict(modelD,phenotypeDataNewWT);
predSlopeNMDA = predict(modelD,phenotypeDataNewNMDA);
predSlopeWT = reshape(predSlopeWT,size(GLU));
predSlopeNMDA = reshape(predSlopeNMDA,size(GLU));
close all
mesh(linspace(0,5000,50)/1000,linspace(0,5000,50)/1000,predSlopeWT);
hold on
mesh(linspace(0,5000,50)/1000,linspace(0,5000,50)/1000,predSlopeNMDA);
xlabel('GLU')
ylabel('GLY')
%%

phenotypeDataNewSWEEP = table;
CON = linspace(0,1000,5);
genoType = 'wt';
growout = '1';
Gly_Conc_um = 1000;
for e = 1:size(CON,2)
    phenotypeDataNewSWEEP{e,'genotype'} = {genoType};
    phenotypeDataNewSWEEP{e,'grow_out'} = {growout};
    phenotypeDataNewSWEEP{e,'Gly_Conc_um'} = Gly_Conc_um;
    phenotypeDataNewSWEEP{e,'Glu_Conc_um'} = CON(e);
end
predSlopeSWEEP = predict(model{toUSE},phenotypeDataNewSWEEP);
predSlopeSWEEPD = predict(modelD,phenotypeDataNewSWEEP);

dDELTA = PCA_BKPROJ(predSlopeSWEEPD,dE,dU);
TIME = linspace(0,numel(dE),numel(dE));
phP = predSlopeSWEEP*TIME + dDELTA;
close all
plot(phP');
figure;
plot(dDELTA');
%%
close all
ypred = predict(complexModel);
plot(ypred,pC(:,1),'.')

%%
model1 = fitlme(phenotypeData,...
'PC_1 ~ 1  + Glu_Conc_um + Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'Distribution','Normal','Link','identity','FitMethod','MPL');

model2 = fitlme(phenotypeData,...
'PC_1 ~ 1  + Glu_Conc_um: Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'Distribution','Normal','Link','identity','FitMethod','MPL');
%%
%%

glme = fitglme(phenotypeData,...
'PC_1 ~ 1  + Glu_Conc_um + Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'Distribution','Normal','Link','identity','FitMethod','MPL');


glme = fitglme(phenotypeData,...
'PC_1 ~ 1  + Glu_Conc_um + Gly_Conc_um + (1|genotype) + (1|grow_out)',...
'Distribution','Normal','Link','identity','FitMethod','MPL');

%%
glme = fitglme(mfr,...
'defects ~ 1 + newprocess + time_dev + temp_dev + supplier + (1|factory)',...
'Distribution','Poisson','Link','log','FitMethod','Laplace',...
'DummyVarCoding','effects')
%% init 
phenotype = bsxfun(@minus,phenotype,phenotype(:,1));
%% compare min to 100um - 100 um GLU ONLY
sel{1} = find(strcmp(amount,'0mm'));
sel{2} = find(strcmp(amount,'100um'));
cnt = 1;
clear U S LEG
%close all
figure
TRL = {'glu' 'gly' 'min'};
for e = 1:numel(sel)
    g = genotype(sel{e});
    UQ = unique(g);
    tmpP = phenotype(sel{e},:);
    tmpA = amount(sel{e});
    tmpT = TYPE(sel{e},:);
    
    for f = 1:numel(UQ)
        sIDX = strcmp(g,UQ{f});
        subP{cnt} = tmpP(sIDX,:);
        U(cnt,:) = mean(subP{cnt},1);
        S(cnt,:) = std(subP{cnt},1)*size(subP{cnt},1).^-.5;
        tempy = find(tmpT(1,:));
       
        LEG{cnt} = [tmpA{1} '--' UQ{f} '--'  [TRL{tempy}]];
        cnt = cnt + 1;
    end
    
end


errorbar(U',S');
legend(LEG)

%% compare min to 100um  - GLY-1mm
sel{1} = find(strcmp(amount,'0mm'));
sel{2} = find(strcmp(amount,'1mm') & all(TYPE == repmat([0 1 0],[size(TYPE,1) 1]),2)');
cnt = 1;
clear U S LEG
figure;
TRL = {'glu' 'gly' 'min'};
for e = 1:numel(sel)
    g = genotype(sel{e});
    UQ = unique(g);
    tmpP = phenotype(sel{e},:);
    tmpA = amount(sel{e});
    tmpT = TYPE(sel{e},:);
    
    for f = 1:numel(UQ)
        sIDX = strcmp(g,UQ{f});
        subP{cnt} = tmpP(sIDX,:);
        U(cnt,:) = mean(subP{cnt},1);
        S(cnt,:) = std(subP{cnt},1)*size(subP{cnt},1).^-.5;
        tempy = find(tmpT(1,:));
       
        LEG{cnt} = [tmpA{1} '--' UQ{f} '--'  [TRL{tempy}]];
        cnt = cnt + 1;
    end
    
end




errorbar(U',S');
legend(LEG)
%% compare min to 100um - 
sel{1} = find(strcmp(amount,'0mm'));
sel{2} = find(strcmp(amount,'1mm') & all(TYPE == repmat([1 1 0],[size(TYPE,1) 1]),2)');
cnt = 1;
clear U S LEG
figure;
TRL = {'glu' 'gly' 'min'};
for e = 1:numel(sel)
    g = genotype(sel{e});
    UQ = unique(g);
    tmpP = phenotype(sel{e},:);
    tmpA = amount(sel{e});
    tmpT = TYPE(sel{e},:);
    
    for f = 1:numel(UQ)
        sIDX = strcmp(g,UQ{f});
        subP{cnt} = tmpP(sIDX,:);
        U(cnt,:) = mean(subP{cnt},1);
        S(cnt,:) = std(subP{cnt},1)*size(subP{cnt},1).^-.5;
        tempy = find(tmpT(1,:));
       
        LEG{cnt} = [tmpA{1} '--' UQ{f} '--'  [TRL{tempy}]];
        cnt = cnt + 1;
    end
    
end


errorbar(U',S');
legend(LEG)
%% compare min to 100um  - GLU-GLY-1mm
sel{1} = find(strcmp(amount,'0mm'));
sel{2} = find(strcmp(amount,'1mm') & all(TYPE == repmat([1 1 0],[size(TYPE,1) 1]),2)');
cnt = 1;
clear U S LEG
figure;
TRL = {'glu' 'gly' 'min'};
for e = 1:numel(sel)
    g = genotype(sel{e});
    UQ = unique(g);
    tmpP = phenotype(sel{e},:);
    tmpA = amount(sel{e});
    tmpT = TYPE(sel{e},:);
    
    for f = 1:numel(UQ)
        sIDX = strcmp(g,UQ{f});
        subP{cnt} = tmpP(sIDX,:);
        U(cnt,:) = mean(subP{cnt},1);
        S(cnt,:) = std(subP{cnt},1)*size(subP{cnt},1).^-.5;
        tempy = find(tmpT(1,:));
       
        LEG{cnt} = [tmpA{1} '--' UQ{f} '--'  [TRL{tempy}]];
        cnt = cnt + 1;
    end
    
end


errorbar(U',S');
legend(LEG)
%% compare min to 100um - GLU -1mm
sel{1} = find(strcmp(amount,'0mm'));
sel{2} = find(strcmp(amount,'1mm') & all(TYPE == repmat([1 0 0],[size(TYPE,1) 1]),2)');
cnt = 1;
clear U S LEG
figure;
TRL = {'glu' 'gly' 'min'};
for e = 1:numel(sel)
    g = genotype(sel{e});
    UQ = unique(g);
    tmpP = phenotype(sel{e},:);
    tmpA = amount(sel{e});
    tmpT = TYPE(sel{e},:);
    %{
    tmpP = imfilter(tmpP,fspecial('average',[1 20]),'replicate');
    tmpP = gradient(tmpP);
        %}
    for f = 1:numel(UQ)
        sIDX = strcmp(g,UQ{f});
        subP{cnt} = tmpP(sIDX,:);
        U(cnt,:) = mean(subP{cnt},1);
        S(cnt,:) = std(subP{cnt},1)*size(subP{cnt},1).^-.5;
        tempy = find(tmpT(1,:));
       
        LEG{cnt} = [tmpA{1} '--' UQ{f} '--'  [TRL{tempy}]];
        cnt = cnt + 1;
    end
    
end


errorbar(U',S');
legend(LEG)
%% bar graphs below
GR = mean(diff(phenotype,1,2),2);
%% compare min to 100um - 100 um GLU ONLY
sel{1} = find(strcmp(amount,'0mm'));
sel{2} = find(strcmp(amount,'100um'));
cnt = 1;
clear U S LEG
%close all
figure
hold on
LEG = {};
TRL = {'glu' 'gly' 'min'};
for e = 1:numel(sel)
    g = genotype(sel{e});
    UQ = unique(g);
    tmpP = GR(sel{e},:);
    tmpA = amount(sel{e});
    tmpT = TYPE(sel{e},:);
    
    for f = 1:numel(UQ)
        sIDX = strcmp(g,UQ{f});
        subP{cnt} = tmpP(sIDX,:);
        U(cnt,:) = mean(subP{cnt},1);
        S(cnt,:) = std(subP{cnt},1)*size(subP{cnt},1).^-.5;
        tempy = find(tmpT(1,:));
        LEG{cnt} = [tmpA{1} '--' UQ{f} '--'  [TRL{tempy}]];
        LEG2{cnt+1} = ['Group ' num2str(cnt)];
        cnt = cnt + 1;
    end
    
end

bar(U');
legend(LEG)
errorb(U',S')
set(gca,'XTickLabel',LEG2)

