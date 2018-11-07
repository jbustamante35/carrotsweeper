%%
%FilePath = '/mnt/spaldingimages/nate/whole_Candace_Gravitropism_arabidopsis/';
%FilePath = '/mnt/scratch1/myLinks/arabidopsisGravitropism/baseLine/';
FilePath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/toProcess/';
FilePath = '/mnt/piranhapuffer2/Hannah/Maize QTL Gravi/rawData/';
%FilePath = '/mnt/scratch1/myLinks/arabidopsisGravitropism/candace_qtl/';
%FilePath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/toProcess/straight/RafaelStraight14Aug13/lip5-1_Cam2/';
%FilePath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/toProcess/12 hours Gravi/RafaelGravi5Sep13/vps2_2_Cam5/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% sort sets
for e = 1:numel(SET)
    NAME = [];
    for i = 1:numel(SET{e})
        [pth,nm,ext] = fileparts(SET{e}{i});
        NAME(i) = str2num(nm);
    end
    [J sidx] = sort(NAME);
    SET{e} = SET{e}(sidx);
    e
end
%% try again
oPath = '/mnt/scratch5/arabidopsis_model_13.09.09/';
oPath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/new_return_13.10.04/';

mkdir(oPath);
for e = 1:numel(SET)
    tm = clock;
    isolateRoots_overStack(SET{e},oPath,1,3,20,20,20,1);
    etm = etime(clock,tm);
    etm
end
%%
oPath = '/mnt/scratch5/arabidopsis_model_13.09.09/';
oPath = '/mnt/piranhapuffer2/Otegui/Rafael/phytoMorph/morphoMetrics/raw_data/new_return_13.10.04/';
oPath = '/mnt/piranhapuffer2/Hannah/Maize QTL Gravi/return/';
mkdir(oPath);
DEPTH = 3;
for e = 1:numel(SET)
    % check for 
    csvOutPath = [oPath 'csv/'];
    mkdir(csvOutPath);
    [pth,nm,ext] = fileparts(SET{e}{1});
    fidx = strfind(pth,filesep);
    pth = pth(fidx(end-DEPTH)+1:end);
    pth = strrep(pth,filesep,'----');
    cFile = [csvOutPath pth '--angle.csv'];
    tm = clock;
    if ~exist(cFile)
        isolateRoots_overStack(SET{e},oPath,1,DEPTH,20,20,20,1);
        etm = etime(clock,tm);
        fprintf(['Stack time is:' num2str(etm) '\n']);
    end
end