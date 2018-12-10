FilePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/return/seedlingData/output/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% treatment scan
key = 'Treatment';
keyID = 'PictureDay';
pKey = 'Phenotype';
clear S
S.hello = 0;
for e = 1:numel(FileList)
    [v] = keyLookup(FileList{e},key);
    id = str2num(keyLookup(FileList{e},keyID));
    pid = keyLookup(FileList{e},pKey);
    v = strrep(v,' ','');
    v = strrep(v,'7','');
    v = strrep(v,'11','');
    if ~isempty(v)
        %if strcmp(pid,'PlantHeight')
        if strcmp(pid,'DigitalBioMass')
            try
                data = mean(csvread(FileList{e}));
                if numel(data) == 1
                    if ~isfield(S,(v))
                        S.(v) = {};
                        S.([v '_n']) = {};
                    end
                    if numel(S.(v)) < id
                        S.(v){id} = [];
                        S.([v '_n']){id} = {};
                    end
                    S.(v){id} = [S.(v){id};data];
                    S.([v '_n']){id}{end+1} =  FileList{e};
                    %waitforbuttonpress
                end
            catch
            end
        end
    end
    e
    numel(FileList)
end
%% treatement scan part two
close all
u = [];
s = [];
for e = 1:21
    u(2,e) = nanmean(S.Control{e});
    u(1,e) = nanmean(S.Cstressday{e});
    s(2,e) = nanstd(S.Control{e})*numel(S.Control{e})^-.5;
    s(1,e) = nanstd(S.Cstressday{e})*numel(S.Cstressday{e})^-.5;
end
errorbar(u',s');
%% find unique types
key = {};
for e = 1:numel(FileList)
    try
        [p,nm,ext] = fileparts(FileList{e});
        n = strfind(nm,'}');
        key{e} = nm(1:n(7));
    catch ME
        if ~isempty(strfind(FileList{e},'.txt'))
            FileList{e}
        end
        [p,nm,ext] = fileparts(FileList{e});
        ext
        key{e} = 'NO KEY';
    end
end
UQ = unique(key);
%% get the unique plots
for e = 1:numel(UQ)
    plotID{e} = keyLookup(UQ{e},'Plot');
end
%% organize by data insert into local "database"
S = [];
for e = 1:numel(FileList)
    e
    numel(FileList)
    try
        plotV = ['plot_' keyLookup(FileList{e},'Plot')];
        expV = keyLookup(FileList{e},'Experiment');
        genoTypeV = keyLookup(FileList{e},'Genotype');
        phenoTypeV = keyLookup(FileList{e},'Phenotype');
        pictureDayV = str2num(keyLookup(FileList{e},'PictureDay'));
        plantNumberV = keyLookup(FileList{e},'PlantNumber');
        if strcmp(plantNumberV,'All')
            plantNumberV = 1:3;
        else
            plantNumberV = str2num(plantNumberV);
        end
        S.(plotV).genoType = genoTypeV;
        S.(plotV).expType = expV;
        if ~strcmp(phenoTypeV,'Curvature')
            S.(plotV).(phenoTypeV).data(pictureDayV,plantNumberV) = csvread(FileList{e});
            S.(plotV).(phenoTypeV).imageFile{pictureDayV}= FileList{e};
        end
    catch
    end
end
%% plot dynamics of Stem-diameter and digital biomass
f = fields(S);
close all
M1 = [];
M2 = [];
M3 = [];
G = {};
T = {};
disp = 0;
for e = 1:numel(f)
    try
        cnt = sum(all(S.(f{e}).StemDiameter.data~=0,2));
        if cnt > 10
            %for k = 1:size(S.(f{e}).StemDiameter.data,2)
                tmp1 = S.(f{e}).StemDiameter.data;
                tmp2 = S.(f{e}).DigitalBioMass.data;
                tmp3 = S.(f{e}).PlantHeight.data;
                tmpG = keyLookup(S.(f{e}).StemDiameter.imageFile{end},'Genotype');
                tmpT = keyLookup(S.(f{e}).StemDiameter.imageFile{end},'Treatment');
                tmp1 = curveFix(tmp1);
                tmp2 = curveFix(tmp2);
                tmp3 = curveFix(tmp3);
                ridx = find(all(tmp1==0,1) | all(tmp2==0,1));
                tmp1(:,ridx) = [];
                tmp2(:,ridx) = [];
                tmp3(:,ridx) = [];
                tmp1 = mean(tmp1,2);
                tmp2 = mean(tmp2,2);
                tmp3 = mean(tmp3,2);
                %tmp1 = imfilter(tmp1,fspecial('average',[3 1]),'replicate');
                %tmp2 = imfilter(tmp2,fspecial('average',[3 1]),'replicate');
                if numel(tmp1) == 22 && numel(tmp2) && numel(tmp3)
                    tmp1(end) = [];
                    tmp2(end) = [];
                    tmp3(end) = [];
                end
                    
                if numel(tmp1) == 21 && numel(tmp2) == 21 && numel(tmp3) == 21
                    M1 = [M1 tmp1];
                    M2 = [M2 tmp2];
                    M3 = [M3 tmp3];
                    G{end+1} = tmpG;
                    T{end+1} = tmpT;
                end
                %{
                rm = find(tmp1==0 | tmp2==0);
                tmp1(rm) = [];
                tmp2(rm) = [];
                %}
                if disp
                    plot(tmp2,tmp1);
                    %hold all
                    drawnow
                    waitforbuttonpress
                end
            %end
        end
        e
        numel(f)
    catch
    end
end
%%
close all
rm = find(any(isnan(M1),1) | any(isinf(M1),1) | any(isnan(M2),1) | any(isinf(M2),1) | any(isnan(M3),1) | any(isinf(M3),1));
M1(:,rm) = [];
M2(:,rm) = [];
M3(:,rm) = [];
T(rm) = [];
G(rm) = [];
LEG = {'B73 Control','B73 Cold','Mo17 Control','Mo17 Cold'}
figure;
plot1(M1,M2,G,'B73',T,'Control','b');
plot1(M1,M2,G,'B73',T,'7C stress day 11','c');
plot1(M1,M2,G,'Mo17',T,'Control','r');
plot1(M1,M2,G,'Mo17',T,'7C stress day 11','m');
legend(LEG)
% make
figure
plotX(M1,G,'B73',T,'Control','b');
plotX(M1,G,'B73',T,'7C stress day 11','c');
plotX(M1,G,'Mo17',T,'Control','r');
plotX(M1,G,'Mo17',T,'7C stress day 11','m');
legend(LEG)
title('Stem Diameter')

figure
plotX(M2,G,'B73',T,'Control','b');
plotX(M2,G,'B73',T,'7C stress day 11','c');
plotX(M2,G,'Mo17',T,'Control','r');
plotX(M2,G,'Mo17',T,'7C stress day 11','m');
legend(LEG)
title('Digital BioMass')

figure
plotX(M3,G,'B73',T,'Control','b');
plotX(M3,G,'B73',T,'7C stress day 11','c');
plotX(M3,G,'Mo17',T,'Control','r');
plotX(M3,G,'Mo17',T,'7C stress day 11','m');
legend(LEG)
title('Plant Height')

%% analysis of M1 M2
rm = find(any(isnan(M1),1) | any(isinf(M1),1) | any(isnan(M2),1) | any(isinf(M2),1));
M1(:,rm) = [];
M2(:,rm) = [];
[S1 C1 U1 E1 L1 ERR1 LAM1] = PCA_FIT_FULL(M1',3);
[S2 C2 U2 E2 L2 ERR2 LAM2] = PCA_FIT_FULL(M2',2);

toE = E1;
toC = C1;
toU = U1;
close all

for s = 1:size(toC,2)
    L = linspace(min(toC(:,s)),max(toC(:,s)),10);
    for l = 1:numel(L)
        C0 = mean(toC,1);
        C0(s) = L(l);
        M = PCA_BKPROJ(C0,toE,toU);
        plot(M)
        hold on
    end
    waitforbuttonpress
    close all
end
corr(C1,C2)
%% plot end point of Stem-diameter and digital biomass
f = fields(S);
close all
for e = 1:numel(f)
    try
        cnt = sum(all(S.(f{e}).StemDiameter.data==0,2));
        plot(S.(f{e}).StemDiameter.data(end,:),S.(f{e}).DigitalBioMass.data(end,:),'.');
        hold all
        drawnow
    catch
        
    end
end
%% plot digital bio mass with curve fix
close all
for e = 1:numel(f)
    try 
        d = S.(f{e}).DigitalBioMass.data;
        [nd] = curveFix(d);
        plot(d,'k');
        hold on
        plot(nd);
        hold off
        pl = f{e};
        title(strrep(pl,'_','**'));
        drawnow
        waitforbuttonpress
    catch
        
    end
end
%%
f = fields(S);
close all
for e = 1:numel(f)
    plot(S.(f{e}).DigitalBioMass.data);
    hold on
    drawnow
end
%%
FUN1_1 = [];
FUN1_2 = [];
%%
oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/return/seedlingData/compiledResults/';
mkdir(oPath);
close all
h1 = figure;
h2 = figure;
toSave = 0;
toDisplay = 0;
M1 = [];
M2 = [];
for u = 1:numel(UQ)
    try
        searchString = UQ{u};
        H = [];
        DBM = [];
        STEM = [];
        pid = '';
        fprintf(['Searching']);
        for e = 1:numel(FileList)
            [p,n,ext] = fileparts(FileList{e});
            fidx = strfind(n,searchString);
            %GENO = keyLookup(FileList{e},'Genotype');
            
            if ~isempty(fidx) %& strcmp(GENO,'B73')
                hello = 1;
                
                pid = keyLookup(FileList{e},'Plot');
                didx = strfind(n,'{PlantNumber_All}{Phenotype_PlantHeight}');
                if ~isempty(didx)
                    gidx = strfind(n,'{PictureDay_');
                    hidx = strfind(n(gidx(1):end),'}');
                    kidx = strfind(n(gidx(1):end),'_');
                    data = csvread(FileList{e});
                    idx = str2num(n((kidx(1)+gidx(1)):(hidx(1)+gidx(1)-2)))-3;
                    H(idx,:) = data;
                end
                
                didx = strfind(n,'{PlantNumber_All}{Phenotype_DigitalBioMass}');
                if ~isempty(didx)
                    gidx = strfind(n,'{PictureDay_');
                    hidx = strfind(n(gidx(1):end),'}');
                    kidx = strfind(n(gidx(1):end),'_');
                    data1 = csvread(FileList{e});
                    idx = str2num(n((kidx(1)+gidx(1)):(hidx(1)+gidx(1)-2)))-3;
                   
                    
                    try
                        if numel(data1)==3
                            for k = 1:3
                                tmpFileName = strrep(FileList{e},'All',num2str(k));
                                tmpFileName = strrep(tmpFileName,'DigitalBioMass','StemDiameter');
                                data = csvread(tmpFileName);
                                STEM(idx,k) = data;
                            end

                            DBM(idx,:) = data1;
                            
                        end
                        
                    catch ME
                        
                        ME
                    end
                    
                end
                
               
                fprintf(['.']);
                if mod(e,30) == 0
                     fprintf(['\n']);
                end 
            end
        end
        
        if ~isempty(DBM) & ~isempty(STEM)
            %{
            FUN1_1 = cat(3,FUN1_1,DBM);
            FUN1_2 = cat(3,FUN1_2,STEM);
            %}
            
            tmp = DBM(end,:);
            M1 = [M1;tmp(:)];
            tmp = STEM(end,:);
            M2 = [M2;tmp(:)];
            %{
            
            idx = find(all(DBM==0,2) & all(STEM==2,2));
            DBM(idx,:) = [];
            STEM(idx,:) = [];
            if size(STEM,1) > 8
                for p = 1:3
                    plot(DBM(:,p),STEM(:,p));
                    drawnow
                    hold on
                end
            end
            %}
            
        end
        
        M1 = [M1;DBM(:)];
        M2 = [M2;STEM(:)];
        
        
        plot(M1(:),M2(:),'.')
        drawnow
        if toDisplay
            figure(h1);
            plot(H)
            oFile = [oPath searchString '{PictureDay_1-22}{PlantNumber_All}{Phenotype_PlantHeight}' '.tif'];
            if toSave;saveas(gca,oFile);end

           

            figure(h2)
            plot(DBM);
            title(num2str(u));
            oFile = [oPath searchString '{PictureDay_1-22}{PlantNumber_All}{Phenotype_DigitalBioMass}' '.tif'];
            if toSave;saveas(gca,oFile);end
            title(pid);

            drawnow
            u
            numel(UQ)




            %waitforbuttonpress
            if toSave
            oFile = [oPath searchString '{PictureDay_1-22}{PlantNumber_All}{Phenotype_PlantHeight}' '.csv'];
            csvwrite(oFile,H);
            oFile = [oPath searchString '{PictureDay_1-22}{PlantNumber_All}{Phenotype_DigitalBioMass}' '.csv'];
            csvwrite(oFile,DBM);
            end
        end
    catch ME
        ME
        FileList{e}
        e
        break
    end
end
