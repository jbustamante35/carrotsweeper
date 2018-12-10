%% Jeff and Nate are using this file as of Dec 03 2015 and of Jan 05 2015
%% 1) scan for tip angle data
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/';
%FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/fullExtract_14.10.14/';
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/fullExtract_12.03.15/';
FileList = {};
FileExt = {'csv'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% 2.X) look for those which have logan spool
ta_idx = [];
le_idx = [];
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'full'))
        if ~isempty(strfind(FileList{e},'angle'))
            ta_idx = [ta_idx ; e];
        end
        if ~isempty(strfind(FileList{e},'length'))
            le_idx = [le_idx ; e];
        end
    
    end
end
AngleFile = FileList(ta_idx);
LengthFile = FileList(le_idx);
%% 2.X) again if only have angles - therefore assign FileList to AngleFile
AngleFile = FileList;
%% 2) look for those which have NOT logan spool
ta_idx = [];
le_idx = [];
for e = 1:numel(FileList)
    if isempty(strfind(FileList{e},'logan')) & isempty(strfind(FileList{e},'length')) & isempty(strfind(FileList{e},'rate'))  
        ta_idx = [ta_idx ; e];
    end
end
AngleFile = FileList(ta_idx);
%% 2.YES) set AngleFile to FileList
AngleFile = FileList;
%% 3) load tip angle
TA = [];
nmTA = [];
rmT = [];

graviPlateName = {};
graviWellName = {};
graviFileName = {};

NOT_READ = {};

for e = 1:numel(AngleFile)
    try 
        tmp = csvread(AngleFile{e});
        if size(tmp,1) == 61
            tmp = tmp';
        end
        
        
        [p,n,ex] = fileparts(AngleFile{e});
        
        
        n = lower(n);
        if strcmp(n(end-3:end),'data')
            n = n(1:end-5);
        end
        
        if ~isempty(strfind(n,'spaldingdata'))
            fidx = strfind(n,'_s');
            n = n((fidx(end)+1):end);
        end
        
        
        
        n = [n '_'];
        fidx = strfind(n,'_');
        tmpPlateName = n((fidx(1)+1):(fidx(2)-1));
        
        fprintf([AngleFile{e} '-->' tmpPlateName '\n']);
        
        
        wellName = n((fidx(2)):end);
        fidx = strfind(wellName,'_');
        wn = {};
        for i = 1:(numel(fidx)-1)
            wn{i} = wellName((fidx(i)+1):(fidx(i+1)-1));
        end


        if numel(wn) == size(tmp,1) & size(tmp,2) == 61;
            
            for i = 1:numel(wn)
                graviPlateName{end+1} = tmpPlateName;
                graviWellName{end+1} = wn{i};
                graviFileName{end+1} = AngleFile{e};
            end
            TA = [TA;tmp];
            nmTA = [nmTA;e*ones(size(tmp,1),1)];
        else
            NOT_READ{end+1} = AngleFile{e};
        end
        TAfn{e} = n;
        
        e
        numel(AngleFile)
   catch ME
       e
       AngleFile{e}
       rmT = [rmT;e];
       
   end
end
%% 4) get unique data
[TA,uidx] = unique(TA,'rows');
graviPlateName = graviPlateName(uidx);
graviWellName = graviWellName(uidx);
graviFileName = graviFileName(uidx);
nmTA = nmTA(uidx);
%% *) search for data files
for e = 1:numel(AngleFile)
    [p,n,ex] = fileparts(AngleFile{e});
    if ~isempty(strfind(n,'22600'))
        AngleFile{e}
        e
    end
end  
%% load length data NOT WORKING YET
LE = [];
nmLE = [];
rmL = [];
for e = 1:numel(LengthFile)
   try
       tmp = csvread(LengthFile{e});
       LE = [LE;tmp];
       nmLE = [nmLE;e*ones(size(tmp,1),1)];
       [p,n,ex] = fileparts(LengthFile{e});
       LEfn{e} = n;
       e
       numel(LengthFile)
   catch ME
       e
       rmL = [rmL;e];
   end
end
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
toRemove = [];
for e = 1:numel(specKey)
    key = specKey{e};
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
specMake(toRM) = [];
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
[kernelVec genoVec popVec plateVec isoVec indexVec femaleVec maleVec] = translateWellNames_forRAW(f);
%% 6.5) average spec together
keys = kernelVec.keySet;
itr = keys.iterator;
import java.util.HashMap;
kernelVecS = HashMap();
kernelVecN = HashMap();
while itr.hasNext()
    tmpKey = itr.next();
    tmpD = kernelVec.get(tmpKey);
    sz = size(tmpD,1);
    fprintf([tmpKey '-->' num2str(sz) '\n']);
    if size(tmpD,2) == 1
        tmpD = tmpD';
    end
    kernelVecS.put(tmpKey,std(tmpD,1,1));
    kernelVecN.put(tmpKey,size(tmpD,1));
    tmpD = mean(tmpD,1);
    kernelVec.put(tmpKey,tmpD);
end
%% JUNK CHECK FOR CALIBRATION
%{
values = popVoc.values();
itr = values.iterator();
while itr.hasNext()
    V{end+1} = 
end
%}
%% 7) START matchup
mSpec = [];
mGravi = [];
mPop = {};
mKey = {};
mGeno = {};
nF = {};
nF_plateName = {};
for e = 1:numel(graviPlateName)
    tmp = graviPlateName{e};
    tmp = makeKey(tmp);
    wn = graviWellName{e};
    wn = upper(wn);
    if numel(wn) < 2
        wn = '00';
    end
    wn = wn(1:2);
    key = [tmp '*_' wn(1:end-1) '_' wn(end)];
    fnd = 0;
    value = kernelVec.get(key);
    popValue = popVec.get(key);
    genoValue = genoVec.get(key);
    
    if ~isempty(value)
        mSpec = [mSpec ;value'];
        mGravi = [mGravi ;TA(e,:)];
        mKey{end+1} = key;
        mPop{end+1} = popValue;
        mGeno{end+1} = genoValue;
        fnd = 1;
    end
    fprintf(['Look up key:' key ':' num2str(fnd) '\n']);
    if fnd == 0
        nF{end+1} = key;
        %nF_plateName{end+1} = graviPlateName{e};
        nF_plateName{end+1} = graviFileName{e};
        
    end
    fprintf(['Done with ' num2str(e) ':' num2str(numel(graviPlateName)) ':' num2str(size(mSpec,1)) '\n']);
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
for i = 6020:numel(vec);%[4380:1:4450]%1:50%[4380:1:4450]%1:numel(total_cross_section_keys)%[4380:1:4450]%[1:50 2000:1:2200 4390:1:4410]% 1:1000%[4390:1:4410]%1:100:numel(total_cross_section_keys)
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
for e = 1:numel(vec)
    fileName = ['/mnt/spaldingdata/nate/JUNK/' num2str(e) '.mat'];
    if exist(fileName)
        load(fileName);
        waitforbuttonpress
        imshow(cIMG,[]);
        hold on
        title(num2str(e))
        rectangle('Position',[BOX_side(1)+1,BOX_side(2)+1,BOX_side(3),BOX_side(4)],'EdgeColor','g');
        rectangle('Position',[d_BOX_front(1)+1,d_BOX_front(2)+1,d_BOX_front(3),d_BOX_front(4)],'EdgeColor','r');
        rectangle('Position',[d_BOX_top(1)+1,d_BOX_top(2)+1,d_BOX_top(3),d_BOX_top(4)],'EdgeColor','b');
        drawnow
    end 
end
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





