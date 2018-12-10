%% scan for files
FilePath = '/mnt/snapper/nate/forRichard/forACA/';
FileList = {};
FileExt = {'nd2'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% label files based on names
wtList = {};
mtList = {};
for e = 1:numel(FileList)
    [p,n,ext] = fileparts(FileList{e});
    n = lower(n);
    if (contains(n,'col') | contains(n,'wt') | contains(n,'ss70')) & ~contains(n,'28') & ~contains(n,'flg22') & ~contains(n,'aca') & ~contains(n,'detached')
        wtList{end+1} = FileList{e};
        %n
    end
    
     if (contains(n,'tl230') | contains(n,'tl231')) & ~contains(n,'28') & ~contains(n,'flg22') & ~contains(n,'detached')
        mtList{end+1} = FileList{e};
        n
    end
end
%% look for "good" data
parfor e = 1:numel(wtList)
    fileName = wtList{e};
    meta = imreadBFmeta(fileName);
    N(e) = meta.nframes;
    %FF(:,:,e) = imreadBF(fileName,1,1,1);
    %imshow(FF(:,:,e),[]);
    %drawnow
end
fidx = find(N < 100);
wtList(fidx) = [];
N(fidx) = [];
%%
oPath = '/mnt/tetra/nate/retWave/wt/';
mkdir(oPath)
for e = 1:numel(wtList)
    wave1(wtList{e},oPath,false);
end

oPath = '/mnt/tetra/nate/retWave/mt/';
mkdir(oPath)
for e = 1:numel(mtList)
    wave1(mtList{e},oPath,false);
end

%% look for "good" data
for e = 1:numel(mtList)
    fileName = mtList{e};
    meta = imreadBFmeta(fileName);
    mN(e) = meta.nframes;
    %mFF(:,:,e) = imreadBF(fileName,1,1,1);
    %imshow(mFF(:,:,e),[70 200]);
    %drawnow
end
fidx = find(mN < 100);
mtList(fidx) = [];
mN(fidx) = [];
%%