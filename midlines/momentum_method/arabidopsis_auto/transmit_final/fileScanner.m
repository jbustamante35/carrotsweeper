%% test with candance
FilePath = '/mnt/spaldingimages/nate/whole_Candace_Gravitropism_arabidopsis/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% loop over the first images for each set
saveStruct.DEPTH = 3;
saveStruct.outPath = '/mnt/scratch5/arabidopsis_model_candace_13.10.15/';
parfor e = 1:numel(SET)    
    out = processImageStack(SET{e},saveStruct,0);
    fprintf(['done with stack@' num2str(e) ':' num2str(numel(SET)) '\n'])
end
%% test with TaN
FilePath = '/mnt/scratch1/myLinks/arabidopsisGravitropism/baseLine/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);