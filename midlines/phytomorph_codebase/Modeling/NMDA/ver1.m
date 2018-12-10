%% to clear all data from matlab workspace
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load from file into a "table" - STEP 1
fileName = '/mnt/snapper/Bessie/NMDAreceptor/NMDA CSV FILES/Length.compiled.final.csv';
T = readtable(fileName);
%% examles for "getting" data from table - STEP NO
% this will list the column names
T.Properties.VariableNames
% this will index into the table row,column
T(3,2)
% index by column name
T.Base(5)
% index many rows
T.Base(3:2:9)
%% load into matrix - STEP 2
phenotype = table2array(T(:,7:end));
%% clean out missing - STEP 3
dP = diff(phenotype,1,2);
mag = 100;
ridx = find(any(abs(dP) > mag*mean(dP(:)),2));
T(ridx,:) = [];
phenotype = table2array(T(:,7:end));
%% normalize data - STEP 4
phenotype = bsxfun(@minus,phenotype,phenotype(:,1));
%% PCA on phenotype - STEP 5
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(phenotype,3);
%% add model values into the table - STEP 6
for e = 1:size(pC,2)
    T{:,['growthPC' num2str(e)]} = pC(:,e);
end
%% obtain the average rate - STEP 7
rate = mean(diff(phenotype,1,2),2);
SIM = rate*(1:size(phenotype,2));
DELTA = phenotype - SIM;
sDELTA = imfilter(DELTA,fspecial('average',[1 21]),'replicate');
[dS dC dU dE dL dERR dLAM] = PCA_FIT_FULL(sDELTA,1);
%% add model values into the table - STEP 8
for e = 1:size(dC,2)
    T{:,['nonlLPC' num2str(e)]} = dC(:,e);
end
%% add average rate - STEP 9
T{:,'averageRate'} = rate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% example model
%model = fitlme(T,'averageRate ~ 1 + Glu:Gly','FitMethod','ML');
model = fitlme(T,'averageRate ~ 1 + Glu:Gly + (1|Growout:Base)','FitMethod','ML');
[~,~,stats] = randomEffects(model)
%% 
modelCompare = compare(model{3},model{1},'nsim',1000,'Options',options);
[~,~,stats] = randomEffects(B);
modelCompare = compare(A,B,'nsim',1000,'Options',options);
%%

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

