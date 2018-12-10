%% Jeff and Nate are using this file as of Feb 4 2016
%% START GRAVI SPOOL

%% STOP GRAVI DATA



%% START SPEC LOAD
%% 1) scan for all csv files for spec data
FilePath = '/mnt/spaldingdata/florida/files/with_time_stamp/';
specFileList = {};
FileExt = {'csv'};
verbose = 1;
specFileList = gdig(FilePath,specFileList,FileExt,verbose);
%% 2) extract extra information and spec from spec files
toTrim = 1679-910+1;
toLoad = numel(specFileList);
%toLoad = 50000;
specData = [];
kernelNumber = [];
specDate = [];
specPacket = [];
specKey = {};
tmSpec = [];

% read the data
tmpData = readtext(specFileList{1});
specData = zeros(toLoad,toTrim);
kernelNumber = zeros(toLoad,1);
specDate = kernelNumber;
specPacket = kernelNumber;
specKey = cell(toLoad,1);
tmSpec = kernelNumber;


parfor e = 1:toLoad
    try
        tm = clock;
        % read the data
        tmpData = readtext(specFileList{e});
        waveLambda = cell2mat(tmpData(2:end,1));
        widx = waveLambda >= 910 & waveLambda <= 1679;
        

        % get meta data from within
        tmpTM = tmpData{1};
        tmpTM = strrep(tmpTM,' ','');
        tmpTM = str2num(tmpTM(2:end-1));
        tmSpec(e) = tmpTM;

        % get meta data from file name
        [pth,nm,ext] = fileparts(specFileList{e});
        sidx = strfind(nm,' ');
        fidx = strfind(nm,'_');
        tmpKN = str2num(nm(1:sidx(1)-1));
        tmpDate = str2num(nm((sidx(1)+1):(fidx(1)-1)));
        tmpPacket = str2num(nm((fidx(1)+1):end));
        tmpSpec = cell2mat(tmpData(widx,2));
        tmpKey = nm((sidx(1)+1):end);

        specData(e,:) = tmpSpec';
        kernelNumber(e) = tmpKN;
        specDate(e) = tmpDate;
        specPacket(e) = tmpPacket;
        specKey{e} = tmpKey;
        dt = etime(clock,tm);

        hrs = (toLoad)*dt/60/60/10;
        fprintf(['loading done in :' num2str(hrs) '\n']);
    
    catch ME
        ME
        e
    end

end
%% 3) remove bad spec loads
for e = 1:numel(specKey)
    toRM(e) = isempty(specKey{e});
end
specData(toRM,:) = [];
kernelNumber(toRM,:) = [];
specDate(toRM,:) = [];
specPacket(toRM,:) = [];
specKey(toRM,:) = [];
tmSpec(toRM,:) = [];
%% STOP SPEC LOAD



%% START WEIGHTS
%% 1) scan for weights files
FilePath = '/mnt/spaldingdata/florida/files/with_time_stamp/';
wtsFileList = {};
FileExt = {'CSV'};
verbose = 1;
wtsFileList = gdig(FilePath,wtsFileList,FileExt,verbose);
%% 2) load wieghts into key,value store
import java.util.Map.*;
import java.util.HashMap;
K2W = HashMap();


toLoad = numel(wtsFileList);

for e = 1:toLoad
    tmpData = readtext(wtsFileList{e});
    [pth,nm,ext] = fileparts(wtsFileList{e});
    nm = strrep(nm,'CornWeights','');
    for i = 1:size(tmpData,1)
        key = [nm '-*-' num2str(tmpData{i,2})];
        value = tmpData{i,3};
        K2W.put(key,value);
    end
    fprintf(['Done with ' num2str(e) ':' num2str(toLoad) '\n']);
end
%% 3) look up of weights for the spec
specWeights = [];
for e = 1:numel(specKey)
    key = [specKey{e} '-*-' num2str(tmSpec(e))];
    value = K2W.get(key);
    specWeights = [specWeights;value];
    fprintf(['Done with ' num2str(e) ':' num2str(numel(specKey)) '\n']);
end
%% STOP WEIGHTS



%% START NIR MASTERLIST
%% 1) load NIR master list
[masterList, result]= readtext('/mnt/spaldingdata/florida/files/nir_master_list.csv');
%% 2) find plate name,genotype for each loaded spec
specGeno = {};
specPlateName = {};
specRepeat = {};
specIsolate = {};
specKernelPosition = {};
specExpCode = {};
specIndex = {};
specFemale = {};
specMale = {};
errorList = {};
specWTS = {};
toRemove = [];
for e = 1:numel(specKey)
    key = specKey{e};
    wkey = [specKey{e} '-*-' num2str(tmSpec(e))];
    fidx = find(strcmp(num2str(key),masterList(:,3)));
    try
        specIndex{end+1} = masterList{fidx,3};
        specGeno{end+1} = masterList{fidx,11}; 
        specPlateName{end+1} = masterList{fidx,7};
        specRepeat{end+1} = masterList(fidx,4);
        specIsolate{end+1} = masterList{fidx,10};
        specKernelPosition{end+1} = masterList{fidx,14};
        specExpCode{end+1} = masterList{fidx,5};
        specFemale{end+1} = masterList{fidx,8};
        specMale{end+1} = masterList{fidx,9};
        specWTS{end+1} = K2W.get(wkey);
    catch ME
        errorList{end+1} = key;
        toRemove = [toRemove e];
    end
    fprintf(['Done with ' num2str(e) ':' num2str(numel(specKey)) '\n']);
end
%% 3) remove bad data sets
toRM = toRemove;
specData(toRM,:) = [];
kernelNumber(toRM,:) = [];
specDate(toRM,:) = [];
specPacket(toRM,:) = [];
specKey(toRM,:) = [];
tmSpec(toRM,:) = [];
specIndex(toRM) = [];
specFemale(toRM) = [];
specMale(toRM) = [];
specWTS(toRM) = [];
if numel(specGeno) ~= numel(specKey)
    specGeno(toRM) = [];
end
if numel(specPlateName) ~= numel(specKey)
    specPlateName(toRM) = [];
end
if numel(specRepeat) ~= numel(specKey)
    specRepeat(toRM) = [];
end
if numel(specIsolate) ~= numel(specKey)
    specIsolate(toRM) = [];
end
if numel(specKernelPosition) ~= numel(specKey)
    specKernelPosition(toRM) = [];
end
if numel(specExpCode) ~= numel(specKey)
    specExpCode(toRM) = [];
end
%% 4) convert specPlateName to strings
for e = 1:numel(specPlateName)
    if ~(ischar(specPlateName{e}))
        specPlateName{e} = num2str(specPlateName{e});
    end
end
%% 5) SIZE CHECK
size(specData)
size(kernelNumber)
size(specDate)
size(specPacket)
size(specKey)
size(tmSpec)
size(specGeno)
size(specPlateName)
size(specIsolate)
size(specKernelPosition)
size(specExpCode)
size(specIndex)
size(specFemale)
size(specMale)
size(specWTS)
%% END NIR MASTERLIST


%% START MATCHING MASTERLIST
%% 6) get key value pair for spec data to match data from gravi
f.specData = specData;
f.POP = specExpCode;
f.position = num2cell(kernelNumber);
f.plateName = specPlateName;
f.genoType = specGeno;
f.specIsolate = specIsolate;
f.specIndex = specIndex;
f.specFemale = specFemale;
f.specMale = specMale;
f.specWTS = specWTS;
[kernelVec genoVec popVec plateVec isoVec indexVec femaleVec maleVec wtsVec] = translateWellNames_forRAW(f);
%% 6.5) average spec together
keys = kernelVec.keySet;
itr = keys.iterator;
import java.util.HashMap;
kernelVecS = HashMap();
kernelVecN = HashMap();
meanStack = [];
stdStack = [];
while itr.hasNext()
    tmpKey = itr.next();
    tmpD = kernelVec.get(tmpKey);
    sz = size(tmpD,1);
    fprintf([tmpKey '-->' num2str(sz) '\n']);
    if size(tmpD,2) == 1
        tmpD = tmpD';
    end
    stmpD = std(tmpD,1,1);
    kernelVecS.put(tmpKey,stmpD);
    kernelVecN.put(tmpKey,size(tmpD,1));
    tmpD = mean(tmpD,1);
    kernelVec.put(tmpKey,tmpD);
    %meanStack = [meanStack;tmpD];
    %stdStack = [stdStack;stmpD];
end
%% check intersection of within and between NOT NEEDED
keySet = kernelVec.keySet;
itr = keySet.iterator();
tmp = kernelVec.get(itr.next());
T1 = zeros(keySet.size(),size(tmp,1));
T2 = zeros(keySet.size(),size(tmp,1));
itr = keySet.iterator();
for e = 1:keySet.size()
    tmpKey = itr.next();
    T1(e,:) = kernelVec.get(tmpKey)';
    T2(e,:) = kernelVecS.get(tmpKey)';
    e
end
[S1 C1 U1 E1 L1 ERR1 LAM1] = PCA_FIT_FULL(T1,3);
[S2 C2 U2 E2 L2 ERR2 LAM2] = PCA_FIT_FULL(T2,3);
[S1,I,S2,F1,F2] = divSpace(E1,E2);
[T1_sub T1_lam] = getSpan(S1,size(S1,2));
[T2_sub T2_lam] = getSpan(S2,size(S1,2));
[I_sub I_lam] = getSpan(I,size(I,2));
%% JUNK CHECK FOR CALIBRATION
%{
values = popVoc.values();
itr = values.iterator();
while itr.hasNext()
    V{end+1} = 
end
%}


%% 7) START matchup of spec and gravi
keySet = kernelVec.keySet;
itr = keySet.iterator();
gravi_data_location = '/mnt/spaldingdata/nate/keyDB_Feb_03_2016/';
cnt = 1;
tcnt = 1;
mKey = {};
mGKey = {};
while (itr.hasNext)
    key = itr.next();
    gravi_key = strrep(key,'*','ASTRIST');
    graviFileName = [gravi_data_location gravi_key '.mat'];
    if exist(graviFileName)
        mKey{cnt} = key;
        mGKey{cnt} = gravi_key;
        cnt = cnt + 1;
    end
    tcnt
    keySet.size()
    tcnt = tcnt + 1;
end
%% remove angle and gr
vecT(:,end-2*61+1:end) = []
%% loader example for gravi quick fun for corr
SD = [];
[vecT SZvec] = loadGraviPhenoTypesByKey([mGKey{1}],50);
vecT = zeros(numel(mGKey),numel(vecT));
for e = 1:numel(mKey)
    tic
    SD(e,:) = mean(kernelVec.get(mKey{e}),1);
    SDe(e,:) = kernelVecS.get(mKey{e});
    tmpG{e} = genoVec.get(mKey{e});
    toc
    tic
    [vecT(e,:) SZvec] = loadGraviPhenoTypesByKey([mGKey{e}],50);
    toc
    e
    numel(mKey)
end
%% quick analysis - decompose
rmidx = any(isnan(vecT),2);
vecT(rmidx,:) = [];
SD(rmidx,:) = [];
SDe(rmidx,:) = [];
[S C U E L ERR LAM] = PCA_FIT_FULL(vecT,5);
[Sx Cx Ux Ex Lx ERRx LAMx] = PCA_FIT_FULL(SD,15);
[Sy Cy Uy Ey Ly ERRy LAMy] = PCA_FIT_FULL(SDe,15);
%% quick analysis - regress
preX = [Cx];
[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(preX,C,5);
%% quick analysis - remove bad data
% make BACKUP of vecT;
%vecTBK = vecT;
%SDBK = SD;
toRM = [1:5];
toRM_N = 60;
RM = [];
for e = 1:numel(toRM)
    [J,sidx] = sort(C(:,toRM(e)));
    RM = [RM sidx(1:toRM_N)];
    RM = [RM sidx(end-toRM_N+1:end)];
end
vecT(RM,:) = [];
SD(RM,:) = [];
SDe(RM,:) = [];
%% quick analysis - sweep
Ypre = [ones(size(preX,1),1) preX]*BETA;
close all
for e = 1:size(Ypre,2)
    figure;
    plot(Ypre(:,e),C(:,e),'.')
end
%% quick display - sweep
close all
uC = mean(C,1);
PER = .8;
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
CL = {'r' 'g' 'b'};
for e = 1:size(C,2)
    
    
    
    tmpC = repmat(uC,[10 1]);
    l = linspace(PER*min(C(:,e)),PER*max(C(:,e)),10);
    
    for k = 1:numel(l)
        tmpC(k,e) = l(k);
    end
    
    
    simM = [];
    fsimM = [];
    %simM = PCA_BKPROJ(C(1:10,:),E,U);
    simM = PCA_BKPROJ(tmpC,E,U);
    for k = 1:numel(l)
        fsimM{k} = graviFold(simM(k,:),SZvec);
    end
    
    for loop = 1:5
        for k = 1:numel(l)
            figure(h1);
            plot(fsimM{k}.kernelContour(:,1,1),-fsimM{k}.kernelContour(:,2,1),CL{1});
            hold on
            plot(fsimM{k}.kernelContour(:,1,30),-fsimM{k}.kernelContour(:,2,30),CL{2});
            plot(fsimM{k}.kernelContour(:,1,61),-fsimM{k}.kernelContour(:,2,61),CL{3});
            plot(fsimM{k}.upperRoot(:,1,1),-fsimM{k}.upperRoot(:,2,1),CL{1});
            plot(fsimM{k}.upperRoot(:,1,30),-fsimM{k}.upperRoot(:,2,30),CL{2});
            plot(fsimM{k}.upperRoot(:,1,61),-fsimM{k}.upperRoot(:,2,61),CL{3});
            plot(fsimM{k}.lowerRoot(:,1,1),-fsimM{k}.lowerRoot(:,2,1),CL{1});
            plot(fsimM{k}.lowerRoot(:,1,30),-fsimM{k}.lowerRoot(:,2,30),CL{2});
            plot(fsimM{k}.lowerRoot(:,1,61),-fsimM{k}.lowerRoot(:,2,61),CL{3});
            plot(fsimM{k}.midline(:,1,1),-fsimM{k}.midline(:,2,1),CL{1});
            plot(fsimM{k}.midline(:,1,30),-fsimM{k}.midline(:,2,30),CL{2});
            plot(fsimM{k}.midline(:,1,61),-fsimM{k}.midline(:,2,61),CL{3});
            axis([-70 250 -100 100]);
            pause(.2);
        end
        hold off
    end
        
    for k = 1:numel(l)  
        
        
         for loop = 1:5
            figure(h1);
            plot(fsimM{k}.kernelContour(:,1,1),-fsimM{k}.kernelContour(:,2,1),CL{1});
            hold on
            plot(fsimM{k}.kernelContour(:,1,30),-fsimM{k}.kernelContour(:,2,30),CL{2});
            plot(fsimM{k}.kernelContour(:,1,61),-fsimM{k}.kernelContour(:,2,61),CL{3});
            plot(fsimM{k}.upperRoot(:,1,1),-fsimM{k}.upperRoot(:,2,1),CL{1});
            plot(fsimM{k}.upperRoot(:,1,30),-fsimM{k}.upperRoot(:,2,30),CL{2});
            plot(fsimM{k}.upperRoot(:,1,61),-fsimM{k}.upperRoot(:,2,61),CL{3});
            plot(fsimM{k}.lowerRoot(:,1,1),-fsimM{k}.lowerRoot(:,2,1),CL{1});
            plot(fsimM{k}.lowerRoot(:,1,30),-fsimM{k}.lowerRoot(:,2,30),CL{2});
            plot(fsimM{k}.lowerRoot(:,1,61),-fsimM{k}.lowerRoot(:,2,61),CL{3});
            plot(fsimM{k}.midline(:,1,1),-fsimM{k}.midline(:,2,1),CL{1});
            plot(fsimM{k}.midline(:,1,30),-fsimM{k}.midline(:,2,30),CL{2});
            plot(fsimM{k}.midline(:,1,61),-fsimM{k}.midline(:,2,61),CL{3});
            axis([-70 200 -100 100]);
        end
        
        for tm = 1:numel(l)
            figure(h2);
            plot(fsimM{tm}.angle,'b');
            hold on
        end
        plot(fsimM{k}.angle,'r');
        hold off
        
        
        
        for tm = 1:numel(l)
            figure(h3);
            gr = mean(gradient(fsimM{tm}.length))*ones(size(fsimM{tm}.length));
            plot(gr,'b');
            hold on
        end
        gr = mean(gradient(fsimM{k}.length))*ones(size(fsimM{tm}.length));
        plot(gr,'r');
        hold off
        
        figure(h4)
        mesh(interp2(-fsimM{k}.kurvature))
        caxis([-.01 .05])
        view([0 90]);
        
        waitforbuttonpress
        figure(h1);
        %hold off
    end
   
    
    
        
end
%% END MATCHING MASTERLIST



%% FEW CORRECTIONS
%% 8) find empty cells and "hand" correct to widiv
% backup BK_POP = mPop
for e = 1:numel(mPop)
    if isempty(mPop{e})
        mPop{e} = 'seedling_phenotyping_widiv';
    end
end
%% 9) get the population counts
UQ = unique(mPop);
TOT = 0;
for u = 1:numel(UQ)
    fidx = strcmp(mPop,UQ{u});
    
    fprintf(['Population count:' UQ{u} ':' num2str(sum(fidx)) '\n']);
    
    %unique(mGeno(fidx))
    
    TOT = TOT + sum(fidx);
end
fprintf(['Population total count:' num2str(TOT) '\n']);
%% 10) LOOK OVER THE NOT FOUND LIST
cnF = {};
cnF_plateName = {};
for e = 1:numel(nF)
    if ~strcmp(nF{e}(1:2),'uw') & ~strcmp(nF{e}(1:2),'ia') & ~strcmp(nF{e}(1:2),'mo') & ~strcmp(nF{e}(1),'p') %& ~strcmp(nF{e}(1:3),'07s')
        if ~isempty(strfind(nF_plateName{e},'Second'))
            cnF{end+1} = nF{e};
            cnF_plateName{end+1} = nF_plateName{e};
        end
        
    end
end
%% FEW CORRECTIONS


%% START SHAPE DATA
% scan for images
FilePath = '/mnt/spaldingdata/KERNELS/imgs/ALL/';
FileList = {};
FileExt = {'tiff'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% add filelist
FilePath2 = '/mnt/spaldingdata/kernelImages/Pictures/';
FileList2 = {};
FileExt = {'tiff'};
verbose = 1;
FileList2 = gdig(FilePath2,FileList2,FileExt,verbose);
%% combine the file lists
for e = 1:numel(FileList2)
    FileList{end+1} = FileList2{e};
    e
    numel(FileList2)
end
%% add filelist
FilePath3 = '/mnt/spaldingdata/KERNELS/imgs/fromKernelImager_December_14_2015/';
FileList3 = {};
FileExt = {'tiff'};
verbose = 1;
FileList3 = gdig(FilePath3,FileList3,FileExt,verbose);
%% combine the file lists
for e = 1:numel(FileList3)
    FileList{end+1} = FileList3{e};
    e
    numel(FileList3)
end
%% backup the file list
FileList_BK = FileList;
%% remove files which have key word color
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'gray'))
        kp(e) = logical(1);
    else
        kp(e) = logical(0);
    end
    e
end
% keep gray images
FileList = FileList(kp);
%% match file names for shape with spec
shape_spec_match = {};
SHAPE_NOT_IN_SPEC = {};
import java.util.HashMap;
topView = HashMap();
sideView = HashMap();
frontView = HashMap();
IND = HashMap();
shape_spec_match_pop = HashMap();
shape_spec_match_gen = HashMap();
for e = 1:numel(FileList)
    fn = FileList{e};
    [pth nm ext] = fileparts(fn);
    fidx = strfind(nm,'_');
    direction{e} = nm(fidx(end)+1:end);
    
    
    tidx = strfind(pth,filesep);
    
    subFolder{e} = pth(tidx(end-1)+1:tidx(end-0)-1);
    
    if numel(fidx) > 1
        
        shapeWellPosition{e} = nm(fidx(end-1)+1:fidx(end)-1);
        shapePlateName{e} = nm(1:fidx(end-1)-1);
        
        
        tmp = shapePlateName{e};
        
        
        
        
        fidx = strfind(tmp,'_');
        if ~isempty(fidx)
            SNIP = tmp(fidx(end)+1:end);
            SNIP = lower(SNIP);
            if (~isempty(strfind(SNIP,'dec')) | ~isempty(strfind(SNIP,'nov')) | ~isempty(strfind(SNIP,'oct')))
                tmp = tmp(1:fidx(end)-1);
            end
            
            if (~isempty(strfind(SNIP,'nor')) | ~isempty(strfind(SNIP,'mut')) | ~isempty(strfind(SNIP,'het')))
                SNIP2 = tmp(fidx(end-1)+1:fidx(end)-1);
                if numel(SNIP2) == 1
                    tmp = strrep(tmp,['_' SNIP2 '_'],['_0' SNIP2 '_']);
                end
            end
            
            
            SNIP = tmp(1:fidx(1)-1);
            SNIP = lower(SNIP);
            
            if (~isempty(strfind(SNIP,'wisn')) & numel(SNIP) == 4)
                tmp(1:4) = [];
                tmp = ['wisn10' tmp];
            end
            
        end
        
        
        
        tmp = makeKey(tmp);
        wn = shapeWellPosition{e};
        
        
        
        if (~isempty(strfind(tmp,'dec')) | ~isempty(strfind(tmp,'nov')) | ~isempty(strfind(tmp,'oct')))
            tmp(end-1:end) = [];
            mtmp = strrep(tmp,'dec','');
            mtmp = strrep(mtmp,'nov','');
            mtmp = strrep(mtmp,'oct','');
            tmp = mtmp;
            tmp(end-1:end) = [];
        end
        
        
        if (~isempty(strfind(tmp,'wisn')))
            if ~strcmp(tmp(end),'s')
                toR = str2num(tmp(end));
                while isempty(toR)
                    tmp
                    tmp(end) = [];
                    toR = str2num(tmp(end));
                    if real(toR) == 0
                        toR = [];
                    end
                end
            end
        end
        
        
        if ~isempty(strfind(tmp,'wisn1015270'))
            tmp = strrep(tmp,'wisn10','wisn08');
        end
        
        
        key = [tmp '*_' wn(1:end-1) '_' wn(end)];
        shapeKey{e} = key;
        
        
        tmp_value = kernelVec.get(key);
        tmp_popinValue = popVec.get(key);
        tmp_genoValue = genoVec.get(key);
        
        
        
        if ~isempty(tmp_value)
            shape_spec_match{end+1} = key;
            % if match then record and save the filename
            if ~isempty(strfind(direction{e},'T'))
                topView.put(key,fn);
            elseif ~isempty(strfind(direction{e},'F')) | ~isempty(strfind(direction{e},'S1'))
                frontView.put(key,fn);
            elseif ~isempty(strfind(direction{e},'S2'))
                sideView.put(key,fn);
            end
            IND.put(key,e);
            shape_spec_match_pop.put(key,tmp_popinValue);
            shape_spec_match_gen.put(key,tmp_genoValue);
        else
            SHAPE_NOT_IN_SPEC{end+1} = key;
        end
        
        
    else
        shapeWellPosition{e} = 'NA';
        shapePlateName{e} = 'NA';
        shapeKey{e} = 'NA';
    end
    e
    numel(FileList)
    %tmp = makeKey(tmp);
end
%% 8) find empty cells and "hand" correct to NC-350-RILS CHECK IF RIGHT
for e = 1:numel(shape_spec_match)
    tmpKey = shape_spec_match{e};
    tmpPop = shape_spec_match_pop.get(tmpKey);
    if isempty(tmpPop)
        shape_spec_match_pop.put(tmpKey,'NC-350_RILs');
    end
end
%% 9) get the population counts
tmpList = shape_spec_match_pop.keySet;
itr = tmpList.iterator();
POPL = {};
KEYL = {};
SUB = {};
while itr.hasNext()
    tmpKey = itr.next();
    POPL{end+1} = shape_spec_match_pop.get(tmpKey);
    KEYL{end+1} = tmpKey;
    tmpIND = IND.get(tmpKey);
    SUB{end+1} = subFolder{tmpIND};
    %POPL{end+1} = itr.next();
end
UQ = unique(POPL);
TOTS = 0;
for u = 1:numel(UQ)
    fidx = strcmp(POPL,UQ{u});
    fprintf(['Population count:' UQ{u} ':' num2str(sum(fidx)) '\n']);
    TOTS = TOTS + sum(fidx);
end
TOTS
%% anaylize image data
close all
disp = 0;
total_cross_section_keys = intersect(shape_spec_match,mKey);
total_cross_section_keys = unique(shape_spec_match);
vec = 1:1:numel(total_cross_section_keys);
ERROR = {};
%vec = fliplr(vec);
for i = 1:numel(vec);%[4380:1:4450]%1:50%[4380:1:4450]%1:numel(total_cross_section_keys)%[4380:1:4450]%[1:50 2000:1:2200 4390:1:4410]% 1:1000%[4390:1:4410]%1:100:numel(total_cross_section_keys)
    try
        e = vec(i);
        tic
        tmpKey = total_cross_section_keys{e};
        fnSIDE = sideView.get(total_cross_section_keys{e});
        fnFRONT = frontView.get(total_cross_section_keys{e});
        fnTOP = topView.get(total_cross_section_keys{e});

        [pth nm1 ext] = fileparts(fnSIDE);
        [pth nm2 ext] = fileparts(fnFRONT);
        [pth nm3 ext] = fileparts(fnTOP);


        fidx = strfind(nm1,'_');
        TYPE1 = nm1(fidx(end)+1:end);

        fidx = strfind(nm2,'_');
        TYPE2 = nm2(fidx(end)+1:end);

        fidx = strfind(nm3,'_');
        TYPE3 = nm3(fidx(end)+1:end);


        [BOX_side d_BOX_front d_BOX_top MASK_side MASK_front MASK_top cIMG flags{i}] = getKERNEL_Data(fnSIDE,fnFRONT,fnTOP);
        if disp
            imshow(cIMG,[]);
            hold on
            rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
            rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
            rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');
            title([TYPE1 '-->' TYPE2 '-->' TYPE3 '-->' num2str(e)]);
            drawnow
            hold off
        end

        save(['/mnt/spaldingdata/nate/JUNK/' num2str(e) '.mat'],'tmpKey','MASK_side','MASK_front','MASK_top','BOX_side','d_BOX_front','d_BOX_top','cIMG'); 
        %{
        SSM{e} = MASK_side;
        SFM{e} = MASK_front;
        STM{e} = MASK_top;
        SSB{e} = BOX_side;
        SFB{e} = d_BOX_front;
        STB{e} = d_BOX_top;
        SCI{e} = cIMG;%[cSIDE cFRONT cTOP];
        %}
        e
        numel(vec)
        toc
    catch
        ERROR{end+1} = i;
    end
    %pause(.1)
end
%% view the shape result from MEORY NOT DISK
close all
for e = 1:numel(vec)
    i = vec(e);
    imshow(SCI{i},[]);
    hold on
    title(num2str(i))
    rectangle('Position',[SSB{i}(1)+1,SSB{i}(2)+1,SSB{i}(3),SSB{i}(4)],'EdgeColor','g');
    rectangle('Position',[SFB{i}(1)+1,SFB{i}(2)+1,SFB{i}(3),SFB{i}(4)],'EdgeColor','r');
    rectangle('Position',[STB{i}(1)+1,STB{i}(2)+1,STB{i}(3),STB{i}(4)],'EdgeColor','b');
    %{
    if sum(SFM{i}(:)) > 1000000
        
        
        fnSIDE = sideView.get(total_cross_section_keys{i});
        fnFRONT = frontView.get(total_cross_section_keys{i});
        fnTOP = topView.get(total_cross_section_keys{i});



        [SSB{i} SFB{i} STM{i} SSM{i} SFM{i} STM{i} SCI{i}] = getKERNEL_Data(fnSIDE,fnFRONT,fnTOP);
        
        
    end
    %}
    hold off
    drawnow
end
%% view the shape results from disk
total_cross_section_keys = intersect(shape_spec_match,mKey);
total_cross_section_keys = unique(shape_spec_match);
vec = 1:1:numel(total_cross_section_keys);
P = [];
Q = [];
N = numel(vec);
%%
FilePath = '/mnt/spaldingdata/nate/JUNK/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% make sheetForJeff
N = 1000;
for e = 1:numel(FileList)
    %fileName = ['/mnt/spaldingdata/nate/JUNK/' num2str(e) '.mat'];
    fileName = FileList{e};
    if exist(fileName)
        load(fileName);
        %waitforbuttonpress
        flag = 0;
        %while( ~flag)
            imshow(cIMG,[]);
            topV = cIMG(:,round(2*size(cIMG,2)/3):end,1:3);
            topV = imresize(topV,size(MASK_top));
            topV = double(topV(:,:,1)).*double(~MASK_top);
            %hold on
            %title([num2str(e) '->' strrep(total_cross_section_keys{e},'_','**')],'FontSize',20)
            rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
            rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
            rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');
            
            %{
            [px py V] = impixel();
            L = size(cIMG,2);
            POS = floor(px/L*3)+1;
            text(px(1),py(1),'Side','Color','w','FontSize',20);
            text(px(2),py(2),'Front','Color','w','FontSize',20);
            button = questdlg('Good Quality?');
            if strcmp(button,'Cancel')
                flag = 0;
            else
                flag = 1;
            end
            Q(e) = strcmp(button,'Yes');
            P = [P;POS'];
            %}
        %end
        drawnow        
        sheetForJeff{e,1} = tmpKey;
        sheetForJeff{e,2} = BOX_side(3);
        sheetForJeff{e,3} = BOX_side(4);
        sheetForJeff{e,4} = d_BOX_front(3);
        sheetForJeff{e,5} = d_BOX_front(4);
        sheetForJeff{e,6} = d_BOX_top(3);
        sheetForJeff{e,7} = d_BOX_top(4);
        sheetForJeff{e,8} = e;
        sheetForJeff{e,9} = fileName;
        sheetForJeff{e,10} = sum(topV(:));
    end 
end
%% look over plate
for e = 1:size(sheetForJeff)
    fidx = strfind(sheetForJeff{e,1},'*');
    KE{e} = sheetForJeff{e,1}(1:fidx(1));
end
%% geneate correction
UQ = unique(KE);
%P = [];
Q = [];
for e = 422:numel(UQ)
    
    
    fidx = find(strcmp(KE,UQ{e}));
    ridx = randperm(numel(fidx));
    
    
    
    %fidx = fidx(ridx(1:min(3,numel(ridx))));
    fidx = fidx(ridx(1));
    
    
    
    for u = 1:numel(fidx)
        fileName = FileList{fidx(u)};
        if exist(fileName)
            load(fileName);
            %waitforbuttonpress
            %flag = 0;
            %while( ~flag)
                imshow(cIMG,[]);
                %hold on
                %title([num2str(e) '->' strrep(total_cross_section_keys{e},'_','**')],'FontSize',20)
                rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
                rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
                rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');

                title(num2str(e));                
                [px py V] = impixel();
                L = size(cIMG,2);
                POS = floor(px/L*3)+1;
                text(px(1),py(1),'Side','Color','w','FontSize',20);
                text(px(2),py(2),'Front','Color','w','FontSize',20);
                %{
                button = questdlg('Good Quality?');
                if strcmp(button,'Cancel')
                    flag = 0;
                else
                    flag = 1;
                end
                Q(e) = strcmp(button,'Yes');
                %}
                P = [P;POS'];
               
            %end
            
            
            
            
            drawnow
            title(UQ{e});
            
            %{
            sheetForJeff{e,1} = tmpKey;
            sheetForJeff{e,2} = BOX_side(3);
            sheetForJeff{e,3} = BOX_side(4);
            sheetForJeff{e,4} = d_BOX_front(3);
            sheetForJeff{e,5} = d_BOX_front(4);
            sheetForJeff{e,6} = d_BOX_top(3);
            sheetForJeff{e,7} = d_BOX_top(4);
            %}
        end 
    end
end
%% find plates with correction signature
close all
toM = [1 2];
toM = [3 3];
sidx = find(~all(bsxfun(@eq,P,toM),2));
sidx = find(all(bsxfun(@eq,P,toM),2));

for e = 1:numel(sidx)
    toMU = UQ(sidx(e));
    
    fidx = find(strcmp(KE,UQ{sidx(e)}));
    P(sidx(e),:)
    for f = 1:numel(fidx)
        fileName = FileList{fidx(f)};
        if exist(fileName)
            load(fileName);
            %waitforbuttonpress
            %flag = 0;
            %while( ~flag)
            imshow(cIMG,[]);
            %hold on
            %title([num2str(e) '->' strrep(total_cross_section_keys{e},'_','**')],'FontSize',20)
            rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
            rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
            rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');

                
            
            
            tidx = strfind(tmpKey,'*');
            wellN = tmpKey(tidx(1)+1:end);
            drawnow
            title([UQ{sidx(e)} '-->' strrep(wellN,'_','')]);
            drawnow
            waitforbuttonpress
            %{
            sheetForJeff{e,1} = tmpKey;
            sheetForJeff{e,2} = BOX_side(3);
            sheetForJeff{e,3} = BOX_side(4);
            sheetForJeff{e,4} = d_BOX_front(3);
            sheetForJeff{e,5} = d_BOX_front(4);
            sheetForJeff{e,6} = d_BOX_top(3);
            sheetForJeff{e,7} = d_BOX_top(4);
            %}
        end
    end
end
%% correct [3 3] signature in plate
toM = [1 2];
toM = [3 3];
sidx = find(~all(bsxfun(@eq,P,toM),2));
sidx = find(all(bsxfun(@eq,P,toM),2));
cP = P;
for e = 1:numel(sidx)
    cP(sidx(e),:) = [1 2];
end
%% write out sheet for jeff
cell2csv('/mnt/spaldingdata/nate/sheetForJeff.csv',sheetForJeff)
%% write out correction data
cell2csv('/mnt/spaldingdata/nate/PVEC.csv',P)
%% list to correct
% 1) empty kernel images
%% make backup of sheetForJeff and fix order of data csheetforJeff
clear csheetForJeff
for k = 1:size(cP,1)
    
    tmp = [cP(k,:) setdiff(1:3,cP(k,:))];
    
    key = UQ{k};
    
    
    
    

    isP = strfind(sheetForJeff(:,1),key);
    cidx = [];
    for i = 1:numel(isP)
        if ~isempty(isP{i})
            cidx = [cidx i];
        end
    end
    
    for j = 1:numel(cidx)
        e = cidx(j);

        
        csheetForJeff{e,1} = sheetForJeff{e,1}; 


        idx1 = [2*tmp(1) 2*tmp(1)+1];

        csheetForJeff{e,2} = sheetForJeff{e,idx1(1)};
        csheetForJeff{e,3} = sheetForJeff{e,idx1(2)};

        idx1 = [2*tmp(2) 2*tmp(2)+1];

        csheetForJeff{e,4} = sheetForJeff{e,idx1(1)};
        csheetForJeff{e,5} = sheetForJeff{e,idx1(2)};    

        idx1 = [2*tmp(3) 2*tmp(3)+1];

        csheetForJeff{e,6} = sheetForJeff{e,idx1(1)};
        csheetForJeff{e,7} = sheetForJeff{e,idx1(2)};    

        csheetForJeff{e,8} = sheetForJeff{e,8};    
        csheetForJeff{e,9} = sheetForJeff{e,9}; 
        
        tmpKey = sheetForJeff{e,1};
        tidx = strfind(tmpKey,'*');
        csheetForJeff{e,10} = tmpKey(1:tidx(1));
        csheetForJeff{e,11} = tmpKey((tidx(1)+1):end);
        
        % background
        csheetForJeff{e,12} = sheetForJeff{e,10};
    end
   
end
%% look for difference in weights between two groups
for e = 1:numel(qidx1)
    try
        d = wtsVec.get(csheetForJeff{qidx1(e),1});
        tmpy1(e) = mean(d);
        %fprintf(['Done with ' num2str(e) ':' num2str(numel(specKey)) '\n']);
    catch
        csheetForJeff{qidx1(e),1}
    end
end
for e = 1:numel(qidx2)
    try
        d = wtsVec.get(csheetForJeff{qidx2(e),1});
        tmpy2(e) = mean(d);
        %fprintf(['Done with ' num2str(e) ':' num2str(numel(specKey)) '\n']);
    catch
        csheetForJeff{qidx2(e),1}
    end
end
close all
[f1,x1] = ksdensity(tmpy1);
[f2,x2] = ksdensity(tmpy2);
plot(x1,f1,'r');
hold on
plot(x2,f2);
(mean(tmpy1))/(mean(tmpy2))
(mean(tmpy1)^(1/3))/(mean(tmpy2)^(1/3))
Rb = (mean(tmpy1.^(1/3)))/(mean(tmpy2.^(1/3)))

Sb = std(tmpy1)^(1/3)/std(tmpy2)^(1/3)
%% explore outliers MASTER PLOT
close all
h1 = cell2mat(csheetForJeff(:,3));
h2 = cell2mat(csheetForJeff(:,5));
h3 = cell2mat(csheetForJeff(:,2));
h4 = cell2mat(csheetForJeff(:,4));
qidx2 = setdiff(1:numel(h1),qidx1)';
qidx3 = find(h1==h2);
qidx2 = setdiff(qidx2,qidx3);
qidx1 = setdiff(qidx1,qidxO1);
qidx2 = setdiff(qidx2,qidxO1);
store = h1(qidx3);

fu1 = intersect(qidx1,subG);
fu2 = setdiff(qidx1,fu1);
% use ratio
%h1 = h1*R1;
%h2 = h2*R2;
%h1 = polyval(p3,h1);
%h2 = polyval(p4,h2);

nm = csheetForJeff(:,1);
fnm = csheetForJeff(:,end-2);

%rmidx = h2 > 600 | h1 > 600;

h2(qidx1) = (h2(qidx1) - b(1))*b(2)^-1;
h3(qidx1) = (h3(qidx1) - b(1))*b(2)^-1;

h2(qidx2) = (h2(qidx2) - b2(1))*b2(2)^-1;
h4(qidx2) = (h4(qidx2) - b2(1))*b2(2)^-1;

% correct ish - via measurements
%h1(qidx2) = h1(qidx2)*R1/R3;
%h2(qidx2) = h2(qidx2)*R1/R3;


newR = .5*(R1+R1N)/(.5*(R3+R3N));
h1(qidx2) = h1(qidx2)*newR;
h2(qidx2) = h2(qidx2)*newR;

h3(qidx2) = h3(qidx2)*newR;
h4(qidx2) = h4(qidx2)*newR;

h1(qidxO1) = -1;
h2(qidxO1) = -1;
h3(qidxO1) = -1;
h4(qidxO1) = -1;



%h1(rmidx) = -1;
%h2(rmidx) = -1;

%
%h1(qidx2) = h1(qidx2)*Rc^-1;
%h2(qidx2) = h2(qidx2)*Rc^-1;
% wrong
%h1(qidx2) = h1(qidx2)*R3/R1;
%h2(qidx2) = h2(qidx2)*R3/R1;

MR = mean(mean([h1(qidx1) h2(qidx1)],2))/mean(mean([h1(qidx2) h2(qidx2)],2))



dd1 = (mean([h1(qidx1) h2(qidx1)],2));
dd2 = (mean([h1(qidx2) h2(qidx2)],2));

figure;
[f1,x1] = ksdensity(dd1);
[f2,x2] = ksdensity(dd2);
plot(x1,f1,'r');
hold on
plot(x2,f2);

figure;

%Rc = (MR.^-1)*Rb;

%h2(qidx1) = (h2(qidx1) - b(1))*(R4/R3)^-1;
%h2(qidx2) = h2(qidx2) - b(1);
%h2(qidx1) = h2(qidx1) - b(1);

%{
h1(qidx1) = h1(qidx1)*R3;
h1(qidx2) = h1(qidx2)*R1;
h2(qidx1) = h2(qidx1)*R4;
h2(qidx2) = h2(qidx2)*R2;
%}
%h2(qidx2) = (h2(qidx2) - b2(1))*b2(2)^-1;

%nm(rmidx) = -1;
%fnm(rmidx) = -1;



%fidx = find(h2 > 400 & h1 < 300);
%close all
plot(h1,h2,'.');

hold on;

plot3(h1(subG),h2(subG),1*ones(numel(subG),1),'r.');
plot3(h1(subG2),h2(subG2),2*ones(numel(subG2),1),'k.');
plot3(h1(subG3),h2(subG3),3*ones(numel(subG3),1),'b.');

%plot(h1(fidx),h2(fidx),'r.')
%plot(h1(qidx1),h2(qidx1),'r.')
%plot(h1(fidx(getString)),h2(fidx(getString)),'r.')
plot([0 1000],[0 1000],'r')

%plot(h1(qidx),h2(qidx),'g.')
%plot(h2(qidx),h1(qidx),'k.')
%{
sel = find(all(bsxfun(@eq,[h1 h2],[432 368]),2));
sel = find(all(bsxfun(@eq,[h1 h2],[374 280]),2));
sel = find(all(bsxfun(@eq,[h1 h2],[268 195]),2));
%sel = 1;
%fileName = fnm{fidx(sel)};
fileName = fnm{sel};
load(fileName);
%waitforbuttonpress
%flag = 0;
%while( ~flag)
figure
imshow(cIMG,[]);

%hold on
%title([num2str(e) '->' strrep(total_cross_section_keys{e},'_','**')],'FontSize',20)
rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');
%}
%% FROM HERE TO HERE 2 for res corrections
%% cut out lasso
figure;
[q1 q2 qidx] = lasso(h1,h2);
%% lasso outs
figure;
[q1 q2 qidxO1] = lasso(h1,h2);
%% view lasso
close all
for e = 1:numel(qidx)
    fileName = fnm{qidx(e)};
    load(fileName);
    %waitforbuttonpress
    %flag = 0;
    %while( ~flag)
    figure
    imshow(cIMG,[]);
    %hold on
    %title([num2str(e) '->' strrep(total_cross_section_keys{e},'_','**')],'FontSize',20)
    
    
    rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
    rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
    rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),BOX_side(4)],'EdgeColor','c');
    rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');
    waitforbuttonpress
    close all
    
end
%% correct lasso for bottom cluster
for e = 1:numel(qidx)
    csheetForJeff{qidx(e),5} = csheetForJeff{qidx(e),3}; 
end
%% correct lasso middle cluster
UQ1 = unique(csheetForJeff(qidx,end-1));
[JUNK j1,j2] = intersect(UQ,UQ1);
for r = 1:numel(j1)
    cP(j1(r),:) = [2 1];
end
%% hand click for resolution correction - scan2
secondImageFile = '/mnt/spaldingdata/Bessie/kernel scans from flatbed/WisconsinDiversityPanel/03May13/WISN10_13757_WISN10_13352_WISN11_22647.tif';
scan2 = imread(secondImageFile);
scan22 = imcrop(scan2);
scan2 = imcrop(scan2);
%% fit lasso
close all
%figure;
%[q1 q2 qidx1] = lasso(h1,h2);
b = robustfit(h1(qidx1),h2(qidx1));
%b = polyfit(h1(qidx1),h2(qidx1),1);
KX = linspace(min(h1(qidx1)),max(h1(qidx1)));
K = polyval(fliplr(b'),KX);
hold on
plot(h1(qidx1),h2(qidx1),'b.')
plot(KX,b(2)^-1*(K-b(1)),'r')
plot([0 1000],[0 1000],'g')
%% fit lasso - group2
close all
figure;
[q1 q2 qidx2s] = lasso(h1,h2);
b2 = robustfit(h1(qidx2s),h2(qidx2s));
%b = polyfit(h1(qidx1),h2(qidx1),1);
KX = linspace(min(h1(qidx2s)),max(h1(qidx2s)));
K = polyval(fliplr(b2'),KX);
figure;
hold on
plot(h1(qidx2s),h2(qidx2s),'b.')
plot(KX,b2(2)^-1*(K-b2(1)),'r')
plot([0 1000],[0 1000],'g')
%% hand click for resolution correction
secondImageFile = '/mnt/spaldingdata/Bessie/kernel scans from flatbed/WisconsinDiversityPanel/03May13/WISN10_13757_WISN10_13352_WISN11_22647.tif';
scan2 = imread(secondImageFile);
scan2 = imcrop(scan2);
%% scan 1
firstImageFile = '/mnt/spaldingdata/Bessie/kernel scans from flatbed/NAM parental lines/07S_2073_05B_07S_2072_04B_07S_2047_01B.tif';
scanR = imread(firstImageFile);
scan12 = imcrop(scanR);
scan1 = imcrop(scanR);
%% click on 
toS = [strrep(lower('07S-2073-05-B'),'-','') '*']
toS = [strrep(lower('07S-2072-04-B'),'-','') '*']
fidx = find(strcmp(csheetForJeff(:,end-1),toS));
nmG = csheetForJeff(fidx,end-2);
nmG1 = csheetForJeff(fidx,1);
for f = 1:numel(nmG)
    fileName = nmG{f};
    load(fileName);
    nmG1(f)
    %waitforbuttonpress
    %flag = 0;
    %while( ~flag)
    figure
    imshow(cIMG,[]);
    %hold on
    %title([num2str(e) '->' strrep(total_cross_section_keys{e},'_','**')],'FontSize',20)
    
    
    rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
    rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
    rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),BOX_side(4)],'EdgeColor','c');
    rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');
    
    figure;
    imshow(scan1)
    waitforbuttonpress
    close all
end
%609-->54254-->07S-2073-05-B*_A_3--NAM_parents
%587-->54255-->07S-2073-05-B*_A_4--NAM_parents
%519-->54257-->07S-2073-05-B*_A_6--NAM_parents
%557-->54261-->07S-2073-05-B*_B_2--NAM_parents
%605-->54262-->07S-2073-05-B*_B_3--NAM_parents
%584-->54265-->07S-2072-04-B*_B_5--NAM_parents
%538-->54269-->07S-2073-05-B*_C_1--NAM_parents
%472-->54272-->07S-2073-05-B*_C_4--NAM_parents
%% look for links
toS = [strrep(lower('WISN10_13352'),'_','') '*']
fidx = find(strcmp(csheetForJeff(:,end-1),toS));
nmG = csheetForJeff(fidx,end-2);
nmG1 = csheetForJeff(fidx,1);
for f = 1:numel(nmG)
    fileName = nmG{f};
    load(fileName);
    nmG1(f)
    %waitforbuttonpress
    %flag = 0;
    %while( ~flag)
    figure
    imshow(cIMG,[]);
    %hold on
    %title([num2str(e) '->' strrep(total_cross_section_keys{e},'_','**')],'FontSize',20)
    
    
    rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
    rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
    rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),BOX_side(4)],'EdgeColor','c');
    rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');
    
    figure;
    imshow(scan2)
    waitforbuttonpress
    close all
end
%% gather clicks for side view
[c13 c23 V] = impixel(scan2);
WID3 = [c13 c23];
dW3 = diff(WID3,1,1);
dW3 = sum(dW3.*dW3,2).^.5;
dW3 = dW3(1:2:end);
%% gather clicks
[c13N c23N V] = impixel(scan22);
WID3N = [c13N c23N];
dW3N = diff(WID3N,1,1);
dW3N = sum(dW3N.*dW3N,2).^.5;
dW3N = dW3N(1:2:end);
%% gather clicks for front view
[c14 c24 V] = impixel(scan2);
WID4 = [c14 c24];
dW4 = diff(WID4,1,1);
dW4 = sum(dW4.*dW4,2).^.5;
dW4 = dW4(1:2:end);
%609-->54254-->07S-2073-05-B*_A_3--NAM_parents
%587-->54255-->07S-2073-05-B*_A_4--NAM_parents
%519-->54257-->07S-2073-05-B*_A_6--NAM_parents
%557-->54261-->07S-2073-05-B*_B_2--NAM_parents
%605-->54262-->07S-2073-05-B*_B_3--NAM_parents
%584-->54265-->07S-2072-04-B*_B_5--NAM_parents
%538-->54269-->07S-2073-05-B*_C_1--NAM_parents
%472-->54272-->07S-2073-05-B*_C_4--NAM_parents
%% gather clicks for side view
[c12 c22 V] = impixel(scan12);
WID2 = [c12 c22];
dW2 = diff(WID2,1,1);
dW2 = sum(dW2.*dW2,2).^.5;
dW2 = dW2(1:2:end);
%% another plate
[c12N c22N V] = impixel(scan1);
WID2N = [c12N c22N];
dW2N = diff(WID2N,1,1);
dW2N = sum(dW2N.*dW2N,2).^.5;
dW2N = dW2N(1:2:end);
%% gather clicks for front view
[c11 c21 V] = impixel(scan1);
WID1 = [c11 c21];
dW1 = diff(WID1,1,1);
dW1 = sum(dW1.*dW1,2).^.5;
dW1 = dW1(1:2:end);
%% make Ratios
toS = [strrep(lower('WISN10_13352'),'_','') '*']
toS0 = [strrep(lower('WISN11_22647'),'_','') '*']
fidx = find(strcmp(csheetForJeff(:,end-1),toS));
fidx0 = find(strcmp(csheetForJeff(:,end-1),toS0));
toS2 = [strrep(lower('07S-2073-05-B'),'-','') '*']
toS3 = [strrep(lower('07S-2072-04-B'),'-','') '*']
fidx2 = find(strcmp(csheetForJeff(:,end-1),toS2));
fidx3 = find(strcmp(csheetForJeff(:,end-1),toS3));
getString2 = [3 4 6 10 11 13 17 20];
getString3 = [2 3 4 5 9 10 11 12 14 17 19 22 26 27 28 29 30 31];
getString = [6 10 11 14];
getString0 = [1 2 6 8 15];
csheetForJeff(fidx(getString),1)
csheetForJeff(fidx2(getString2),1)
csheetForJeff(fidx3(getString3),1)
csheetForJeff(fidx0(getString0),1)
%{
tX1 = cell2mat(csheetForJeff(fidx(getString),2));
tX2 = dW3;
p3 = polyfit(tX1,tX2,1);

tX1 = cell2mat(csheetForJeff(fidx(getString),4))-b(1);
tX2 = dW4;
p4 = polyfit(tX1,tX2,1);
p4(2) = 0;

tX1 = cell2mat(csheetForJeff(fidx2(getString2),2));
tX2 = dW1;
p1 = polyfit(tX1,tX2,1);

tX1 = cell2mat(csheetForJeff(fidx2(getString2),4))-b(1);
tX2 = dW2;
p4 = polyfit(tX1,tX2,1);
p2(2) = 0;
%}


tX1 = cell2mat(csheetForJeff(fidx2(getString2),2));
tX2 = dW2;
p1 = robustfit(tX1,tX2);

tX1 = cell2mat(csheetForJeff(fidx3(getString3),2));
tX2 = dW2N;
p12 = robustfit(tX1,tX2);

tX1 = [cell2mat(csheetForJeff(fidx2(getString2),2));cell2mat(csheetForJeff(fidx3(getString3),2))];
tX2 = [dW2;dW2N];
pT = polyfit(tX1,tX2,1);




R3 = dW3.*(cell2mat(csheetForJeff(fidx(getString),2))-0).^-1;
R3 = mean(R3);

R3N = dW3N.*(cell2mat(csheetForJeff(fidx0(getString0),2))-0).^-1;
R3N = mean(R3N);

R4 = dW4.*(cell2mat(csheetForJeff(fidx(getString),4))-b(1)).^-1;
R4 = mean(R4);


R1 = dW2.*(cell2mat(csheetForJeff(fidx2(getString2),2))-0).^-1;
R1 = mean(R1);
close all
plot(cell2mat(csheetForJeff(fidx3(getString3),2)),dW2N,'.')
hold on
plot(cell2mat(csheetForJeff(fidx2(getString2),2)),dW2,'r.')

R1N = dW2N.*(cell2mat(csheetForJeff(fidx3(getString3),2))-0).^-1;
R1N = mean(R1N);

R2 = dW1.*(cell2mat(csheetForJeff(fidx2(getString2),4))-b(1)).^-1;
R2 = mean(R2);
%% FROM HERE TO HERE 2 for res corrections

%% !! after master plot !! master plot version TOP
d1 = cell2mat(csheetForJeff(:,2));
d2 = cell2mat(csheetForJeff(:,6));
%d1(qidx2) = d1(qidx2)*newR;
bkg = cell2mat(csheetForJeff(:,12));
thresh = 2.2 * 10^8;
%selidx = [find(bkg < thresh); qidx2];
selidx = [find(bkg < thresh)];
selidx = setdiff(1:size(d1,1),selidx);
%close all
figure;
plot(bkg(selidx),'.');
figure
plot(d1(selidx),d2(selidx),'.');

figure;
[q1 q2 TOPidx1] = lasso(d1(selidx),d2(selidx));

%% lasso noise
[q1 q2 NOIidx] = lasso(d1(selidx),d2(selidx));
NOIidx = selidx(NOIidx);
%% look up dates for groups
clear tmpY_1 tmpM_1 tmpD_1
for e = 1:numel(subG)
    key = csheetForJeff{subG(e),1};
    value = indexVec.get(key);
    tmpY_1(e) = str2num(value(1:4));
    tmpM_1(e) = str2num(value(5:6));
    tmpD_1(e) = str2num(value(7:8));
    DN1(e) = datenum(tmpY_1(e),  tmpM_1(e), tmpD_1(e));
end

clear tmpY_2 tmpM_2 tmpD_2
for e = 1:numel(subG2)
    key = csheetForJeff{subG2(e),1};
    value = indexVec.get(key);
    tmpY_2(e) = str2num(value(1:4));
    tmpM_2(e) = str2num(value(5:6));
    tmpD_2(e) = str2num(value(7:8));
    DN2(e) = datenum(tmpY_2(e),  tmpM_2(e), tmpD_2(e));
end

clear tmpY_3 tmpM_3 tmpD_3
for e = 1:numel(subG3)
    key = csheetForJeff{subG3(e),1};
    value = indexVec.get(key);
    tmpY_3(e) = str2num(value(1:4));
    tmpM_3(e) = str2num(value(5:6));
    tmpD_3(e) = str2num(value(7:8));
    DN3(e) = datenum(tmpY_3(e),  tmpM_3(e), tmpD_3(e));
end
%% plot dates
figure
X1 = 1:numel(DN3);
plot(X1,sort(DN3),'r')
hold on
X1 = (numel(DN3)+1):(numel(DN3)+1+numel(DN2));
plot(X1,sort(DN2),'g')
X1 = (numel(DN3)+1):(numel(DN3)+1+numel(DN2));
plot(X1,sort(DN1),'b')

Y1 = datenum(2012,1,1);
Y2 = datenum(2013,1,1);
Y3 = datenum(2014,1,1);
Y4 = datenum(2015,1,1);

plot(Y1*ones(1,numel(DN1)))
plot(Y2*ones(1,numel(DN1)))
plot(Y3*ones(1,numel(DN1)))
plot(Y4*ones(1,numel(DN1)))

%% !!!
subG = selidx(TOPidx1);
subG2 = setdiff(selidx,subG);
subG3 = find(bkg < thresh);
% remove noise
subG = setdiff(subG,NOIidx);
subG2 = setdiff(subG2,NOIidx);
subG3 = setdiff(subG3,NOIidx);


figure;
J = 1;
while ~isempty(J)
    [J iIDX1 iIDX2] = intersect(csheetForJeff(subG,10),csheetForJeff(subG2,10))
    subG2 = [subG2 subG(iIDX1)];
    subG(iIDX1) = [];
end

plot(d1(subG2),d2(subG2),'b.','MarkerSize',2);
hold on
plot(d1(subG),d2(subG),'g.','MarkerSize',2);
plot(d1(subG(iIDX1)),d2(subG(iIDX1)),'r.','MarkerSize',5);
plot(d1(subG2(iIDX2)),d2(subG2(iIDX2)),'m.','MarkerSize',5);



figure;
[JCH iIDX1CH iIDX2CH] = intersect(csheetForJeff(subG,10),csheetForJeff(subG2,10))
plot(d1(subG2),d2(subG2),'b.','MarkerSize',2);
hold on
plot(d1(subG),d2(subG),'g.','MarkerSize',2);
plot(d1(subG(iIDX1)),d2(subG(iIDX1)),'r.','MarkerSize',5);
plot(d1(subG2(iIDX2)),d2(subG2(iIDX2)),'m.','MarkerSize',5);

numel(subG2_0)
numel(subG2_0)/numel(subG2)
numel(subG_0)
numel(subG_0)/numel(subG)


csheetForJeff(subG2(iIDX2),1)
csheetForJeff(subG(iIDX1),1)



figure
subG_0 = [];
fidx = strfind(csheetForJeff(subG,1),'wisn101');
for k = 1:numel(fidx)
    if ~isempty(fidx{k})
        subG_0 = [subG_0;k];
    end
end
plot(d1(subG),d2(subG),'.','MarkerSize',2);
hold on
plot(d1(subG(subG_0)),d2(subG(subG_0)),'ro','MarkerSize',2);


figure
subG2_0 = [];
fidx = strfind(csheetForJeff(subG2,1),'wisn101');
for k = 1:numel(fidx)
    if ~isempty(fidx{k})
        subG2_0 = [subG2_0;k];
    end
end
plot(d1(subG2),d2(subG2),'.','MarkerSize',2);
hold on
plot(d1(subG2(subG2_0)),d2(subG2(subG2_0)),'ro','MarkerSize',2);






%% MAKE FINAL CSV MASTER FILE
HEADER = {};
SECTION_0 = 'key';
SECTION_1 = 'SPEC_';
SECTION_2 = 'SHAPE_';
SECTION_3 = 'GRAVI_';
HEADER{1} = [SECTION_0];
% for spec
HEADER{end+1} = [SECTION_1 'genotype'];
HEADER{end+1} = [SECTION_1 'plateName'];
HEADER{end+1} = [SECTION_1 'repeat'];
HEADER{end+1} = [SECTION_1 'isolate'];
HEADER{end+1} = [SECTION_1 'position'];
HEADER{end+1} = [SECTION_1 'expCode'];
HEADER{end+1} = [SECTION_1 'index'];
HEADER{end+1} = [SECTION_1 'female'];
HEADER{end+1} = [SECTION_1 'male'];
HEADER{end+1} = [SECTION_1 'repeats'];
for e = 1:770
    HEADER{end+1} = [SECTION_1 'specAverage-' num2str(e)];
end
for e = 1:770
    HEADER{end+1} = [SECTION_1 'specStd-' num2str(e)];
end
% for shape
HEADER{end+1} = [SECTION_2 'sideFileName'];
HEADER{end+1} = [SECTION_2 'frontFileName'];
HEADER{end+1} = [SECTION_2 'topFileName'];
HEADER{end+1} = [SECTION_2 'plateName'];
HEADER{end+1} = [SECTION_2 'wellPosition'];
HEADER{end+1} = [SECTION_2 'extractedKey'];
HEADER{end+1} = [SECTION_2 'side_depth'];
HEADER{end+1} = [SECTION_2 'side_height'];
HEADER{end+1} = [SECTION_2 'front_width'];
HEADER{end+1} = [SECTION_2 'front_height'];
HEADER{end+1} = [SECTION_2 'top_depth'];
HEADER{end+1} = [SECTION_2 'top_width'];
HEADER{end+1} = [SECTION_2 'matDataFile'];
% for gravi
HEADER{end+1} = [SECTION_3 'folderName'];
HEADER{end+1} = [SECTION_3 'wellPosition'];
HEADER{end+1} = [];
%% look at top view background
for e = 1:1:numel(total_cross_section_keys)
    fnTOP = topView.get(total_cross_section_keys{e});
    tmpT(:,:,e) = imread(fnTOP);
    imshow(tmpT(:,:,e),[]);
    drawnow
end
%% analysis on all images from side
clear tmp;
close all
total_cross_section_keys = intersect(shape_spec_match,mKey);
disp = 1;
Pnew = [];
for e = 1:numel(total_cross_section_keys)
    
    
    fn = sideView.get(total_cross_section_keys{e});
    tmp = imread(fn);
    
    
    [pth nm ext] = fileparts(fn);
   
    if backGroundKernel(tmp) > 255/2 
        tmp = imcomplement(tmp);
    end

    if strcmp(ext,'_S2.tiff')
        tmp(end-245:end,:) = [];    
    else
        tmp(end-100:end,:) = [];    
    end
    
    
    E = edge(tmp,'canny');
    E = bwareaopen(E,10);
    [H,T,R] = hough(E,'Theta',-90);
    P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:)))) ;           
    lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
    max_len = 0;
    xy_long = [];
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];               

       %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       %plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       %plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    BOX = [0 0 size(tmp,2) xy_long(3)];



    level = graythresh(tmp);
    btmp = double(tmp)/255 > level;


    tmp = imcrop(tmp,BOX);
    %imshow(tmp,[])

    PER = .2;
    u1 = mean(tmp,1);            
    s1 = std(double(tmp),1,1);
    s1 = bindVec(s1);
    level = graythresh(s1);
    level = level - PER*level;
    idx1 = find(s1 > level);
    ux = mean(idx1);
    MAX1 = max(idx1);
    MIN1 = min(idx1);


    u2 = mean(tmp,2);            
    s2 = std(double(tmp),1,2);
    s2 = bindVec(s2);
    level = graythresh(s2);
    level = level - PER*level;
    idx2 = find(s2 > level);
    ux = mean(idx2);
    MAX2 = max(idx2);
    MIN2 = min(idx2);
    MAX2 = size(tmp,1);

    if disp
        imshow(tmp,[])
        hold on
        rectangle('Position',[MIN1,MIN2,MAX1-MIN1,MAX2-MIN2],'EdgeColor','r');
        [pth nm ext] = fileparts(fn);
        fidx = strfind(nm,'_');
        title([nm(fidx(end)+1:end) '-->' e]);
        drawnow
        hold off
    end
    
    Pnew(e,:) = [MAX1-MIN1,MAX2-MIN2];
   
    e
    numel(total_cross_section_keys)
    
end
%% analysis on all images from front
close all
ext = '_F.tiff';
PnewF = zeros(numel(f.imageTOP),2);
imageList = f.imageTOP;
imageList = imageList(randperm(numel(imageList)));
clear tmp;
PER = .65;
disp = 1;

UQ = unique(f.genoType);
widx = [];
for u = 1:numel(UQ)
    cnt = 1;
    for e = 1:numel(f.genoType)
        
        if strcmp(f.genoType{e},UQ{u}) & cnt < 5
            widx = [widx e];
            cnt = cnt + 1;
        end        
    end
end
%parfor e = 1:numel(f.imageTOP)
for e = 1:numel(imageList)
    
    if backGroundKernel(tmp) > 255/2 
        tmp = imcomplement(tmp);
    end
    
    
    tmp(end-100:end,:) = [];
    E = edge(tmp,'canny');
    E = bwareaopen(E,10);
    [H,T,R] = hough(E,'Theta',-90);
    P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:)))) ;           
    lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
    %imshow(E,[]);
    %hold on
    max_len = 0;
    xy_long = [];
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];               

       %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       %plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       %plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    BOX = [0 0 size(tmp,2) xy_long(3)];
        

    level = graythresh(tmp);
    btmp = double(tmp)/255 > level;

    if isempty(BOX)
        figure;
        [J BOX] = imcrop(tmp);
        close all;
    else
        %figure(d);
    end

    
    % modif crop BOX by hand
    BOX(4) = BOX(4) - 50;
    tmp = imcrop(tmp,BOX);
    
    
    tmp1 = stdfilt(tmp,ones(11));
    %tmp1 = imfilter(tmp1,fspecial('average',5));
    u1 = mean(tmp1,1);            
    s1 = std(double(tmp1),1,1);
    u1 = bindVec(u1);
    s1 = bindVec(s1);
    %s1 = s1.*u1;    
    s1 = u1;
    s1 = bindVec(s1)+.1;
    s1 = log(s1);
    s1 = bindVec(s1);
    level = graythresh(s1);
    level = level - PER*level;
    level = .10;
    binaryL = s1 > level;
    R = regionprops(binaryL,'Area','PixelIdxList');
    [~,sidx] = sort([R.Area]);
    myThresh = zeros(size(binaryL));
    myThresh(R(sidx(end)).PixelIdxList) = 1;
    
    idx1 = find(myThresh);
    
    ux = mean(idx1);
    MAX1 = max(idx1);
    MIN1 = min(idx1);

    u2 = mean(tmp,2);            
    s2 = std(double(tmp),1,2);
    s2 = bindVec(s2);
    level = graythresh(s2);
    level = level - PER*level;
    idx2 = find(s2 > level);
    ux = mean(idx2);
    MAX2 = max(idx2);
    MIN2 = min(idx2);
    MAX2 = size(tmp,1);


    %RECTVALUE = f.shapeDataL((e),[1 3]);
    if disp
        imshow(tmp1)    
        rectangle('Position',[MIN1,MIN2,MAX1-MIN1,MAX2-MIN2],'EdgeColor','r');
        drawnow
    end
    %BOTLEFT = [MIN1 MAX2];
    %UPPERLEFT = [MIN1 MAX2 - RECTVALUE(2)];
    %rectangle('Position',[UPPERLEFT(1),UPPERLEFT(2),RECTVALUE(1),RECTVALUE(2)],'EdgeColor','g');



    %Pnew = [Pnew;[MAX1-MIN1,MAX2-MIN2]];
    %Pold = [Pold;f.shapeDataL((e),[1 3])];
    PnewF(e,:) = [MAX1-MIN1,MAX2-MIN2];   
    %Pold(e,:) = [f.shapeDataL((e),[1 3])];    
    %drawnow;            
    %sI = [sI tmp];
    %e = e + 1;cnt = cnt +1;
    e
    numel(imageList)
end
PnewF(:,2) = PnewF(:,2) + 50;
%% sniff at the background
BK = [];
BKNM = {};
for e = 1:10:numel(total_cross_section_keys)
    fnSIDE = sideView.get(total_cross_section_keys{e});
    fnFRONT = frontView.get(total_cross_section_keys{e});
    fnTOP = topView.get(total_cross_section_keys{e});
    
    [pth nm1 ext] = fileparts(fnSIDE);
    [pth nm2 ext] = fileparts(fnFRONT);
    [pth nm3 ext] = fileparts(fnTOP);
    
    
    
    fidx = strfind(nm1,'_');
    TYPE1 = nm1(fidx(end)+1:end);
    
    fidx = strfind(nm2,'_');
    TYPE2 = nm2(fidx(end)+1:end);
    
    fidx = strfind(nm3,'_');
    TYPE3 = nm3(fidx(end)+1:end);
    
    tmpS = imread(fnSIDE);
    tmpF = imread(fnFRONT);
    tmpT = imread(fnTOP);
    
    BK = [BK mean([tmpS(:,1:SNIP) tmpS(:,end-SNIP:end)],2)];
    BKNM{end+1} = TYPE1;
    
    
    BK = [BK mean([tmpF(:,1:SNIP) tmpF(:,end-SNIP:end)],2)];
    BKNM{end+1} = TYPE2;
    
    BK = [BK mean([tmpT(:,1:SNIP) tmpT(:,end-SNIP:end)],2)];
    BKNM{end+1} = TYPE3;
    e
    numel(total_cross_section_keys)/10
end
%% cluster the groups of backgrounds
[bS bC bU bE bL bERR bLAM] = PCA_FIT_FULL_T(BK,3);
%% END SHAPE



%% junk list the plate names that are not matched
for e = 1:numel(cF_plateName)
    cF_plateName{e}
end
%% JUNK FIND calibration
fidx = find(strcmp(POPL,'calibration'));
KEYL(fidx)
%% check when folder name and file name do not match
toCheck = 'NAM_parents';
toCheck = 'seedling_phenotyping_widiv';

fidx = find(strcmp(POPL,toCheck));

key_from_file =  KEYL(fidx);
key_from_folder =  SUB(fidx);


for e = 1:numel(key_from_file)
    tmpfileKEY = key_from_file{e}(1:end-5);
    tmpfolderKEY = makeKey(key_from_folder{e});
    
    
    if (~isempty(strfind(tmpfolderKEY,'dec')) | ~isempty(strfind(tmpfolderKEY,'nov')) | ~isempty(strfind(tmpfolderKEY,'oct')))
        tmpfolderKEY(end-1:end) = [];
        mtmp = strrep(tmpfolderKEY,'dec','');
        mtmp = strrep(mtmp,'nov','');
        mtmp = strrep(mtmp,'oct','');
        tmpfolderKEY = mtmp;
        tmpfolderKEY(end-1:end) = [];
    end

    
    if ~ strcmp(tmpfileKEY,tmpfolderKEY)
        e
        tmpfileKEY
        tmpfolderKEY
    end
end


%% check unique plant names for a population
toCheck = 'NAM_parents';
toCheck = 'seedling_phenotyping_widiv';
toCheck = 'composition_mutants';
fidx = find(strcmp(POPL,toCheck));
tmpPL_M = KEYL(fidx);

for e = 1:numel(tmpPL_M)
    tmpPL_M{e} = tmpPL_M{e}(1:end-5);
end
UQ_M = unique(tmpPL_M)
for e = 1:numel(UQ_M)
    fidx = strcmp(tmpPL_M,UQ_M{e});
    fprintf([UQ_M{e} '-->' num2str(sum(fidx)) '\n']);
end
numel(UQ_M)
%% check format errors for plates containing WISN
tmpPL = {};
for e = 1:numel(shapePlateName)  
    tmp = shapePlateName{e};
    
    
    fidx = strfind(tmp,'_');
    if ~isempty(fidx)
        SNIP = tmp(fidx(end)+1:end);
        SNIP = lower(SNIP);
        if (~isempty(strfind(SNIP,'dec')) | ~isempty(strfind(SNIP,'nov')) | ~isempty(strfind(SNIP,'oct')))
            tmp = tmp(1:fidx(end)-1);
        end
    end
    
    tmp = makeKey(tmp);
    
    
        
    if (~isempty(strfind(tmp,'dec')) | ~isempty(strfind(tmp,'nov')) | ~isempty(strfind(tmp,'oct')))
        tmp(end-1:end) = [];
        mtmp = strrep(tmp,'dec','');
        mtmp = strrep(mtmp,'nov','');
        mtmp = strrep(mtmp,'oct','');
        tmp = mtmp;
        tmp(end-1:end) = [];
    end
    
    
    if strcmp(tmp,'wisn101343')
        e
    end
    
    
    if ~isempty(strfind(tmp,'wisn'))
        tmpPL{end+1} = tmp;
    end

end
%% measure the tmpPL length
UQ = unique(tmpPL);
for e = 1:numel(UQ)
    len(e) = numel(UQ{e});
end
fidx = len ~= 11;
UQ(fidx)
%% intsrect matched with wrong length plate names
setdiff(UQ(fidx),UQ_M)


%% look into NAM parents
fidx = strcmp(mPop,UQ{1});

%% look for strcontains
for e = 1:numel(FileList)
    %if ~isempty(strfind(mKey{e},'13757'))
    if ~isempty(strfind(FileList{e},'13757'))
        e
        %mKey{e}
    end
end
%% JUNK cut down list
UQ = unique(SHAPE_NOT_IN_SPEC);
SHORT_LIST = {};
for e = 1:numel(UQ)
    if ~strcmp(UQ{e}(1:2),'uw') & ~strcmp(UQ{e}(1:2),'ia') & ~strcmp(UQ{e}(1:4),'seeed') & ~strcmp(UQ{e}(1:3),'kss') & ~strcmp(UQ{e}(1:3),'kru')
        SHORT_LIST{end+1} = UQ{e};
    end
end
%% junk
for e = 1:numel(shapeKey)
    %tmp = makeKey(specPlateName{e});
    tmp = shapeKey{e};
    if ~isempty(strfind(tmp,'23070'))
        e
        tmp
        %num2cell(shapeKey(e))
    end
end
%% JUNK code LOST NC350 - 17 plates no raw data for tip angle
toTest = specExpCode;
Sfidx = find(strcmp(toTest,UQ{2}));
tmp = specGeno(Sfidx)
k1 = unique(specPlateName(Sfidx));
k2 = unique(graviPlateName);
for e = 1:numel(k1)
    k1{e} = makeKey(k1{e});
end
for e = 1:numel(k2)
    k2{e} = makeKey(k2{e});
end
mp = setdiff(k1,k2);
%%
%% JUNK code LOST WIDIV - 17 plates no raw data for tip angle
toTest = specExpCode;
Sfidx = find(strcmp(toTest,UQ{4}));
tmp = specGeno(Sfidx)
k1 = unique(specPlateName(Sfidx));
k2 = unique(graviPlateName);
for e = 1:numel(k1)
    k1{e} = makeKey(k1{e});
end
for e = 1:numel(k2)
    k2{e} = makeKey(k2{e});
end
mp = setdiff(k1,k2);
%% JUNK OCDEE
toTest = specExpCode;
Sfidx = find(strcmp(mPop,UQ{2}));
tmp = mKey(Sfidx);
mGeno(fidx);
for e = 1:numel(tmp)
    fidx = strfind(tmp{e},'*');
    tmp{e} = tmp{e}(1:fidx(1)-1);
end
UQ2 = unique(tmUp);
for u = 1:1%numel(UQ2)
    l = sum(strcmp(UQ2{u},tmp));
    kidx = find(strcmp(UQ2{u},tmp));
    mKey(Sfidx(kidx))
    fprintf([UQ2{u} '-->' num2str(l) '\n']);
end
%% temp code
kidx = [];
for e = 1:numel(specPlateName)
    if ~isempty(specPlateName{e}) 
        if isempty(strfind(specPlateName{e},'Cali')) & strcmp(lower(specPlateName{e}(1:2)),'07')
            kidx = [kidx e];
        end
    end
end
%%
kidx2 = [];
for e = 1:numel(nF)
    if ~isempty(nF{e}) 
        if isempty(strfind(nF{e},'Cali')) & strcmp(lower(nF{e}(1:2)),'07')
            kidx2 = [kidx2 e];
        end
    end
end
%% 
%specPlateName(kidx(1:200))
nF(kidx2)





