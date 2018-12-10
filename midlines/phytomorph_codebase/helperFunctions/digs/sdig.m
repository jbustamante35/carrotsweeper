function [FileList] = sdig(FilePath,FileList,FileExt,verbose)
    %%% get file list
    tmpList = gdig(FilePath,{},FileExt,verbose);
    [FileList] = orderFrom_gdig(tmpList,FileList);
end

%{
    FilePath = '/mnt/spaldingimages/';
    FileList = {};
    FileExt = {'tiff','TIF'};
    verbose = 1;
    SET = sdig(FilePath,FileList,FileExt,verbose);
%}