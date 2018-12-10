function [FileList] = sdig(FilePath,FileList,FileExt,verbose)
    
    %%% get file list
    tmpList = gdig(FilePath,{},FileExt,verbose);
    pth = {};
    %%% sep into sets
    for i = 1:numel(tmpList)
        [pth{i} nm{i}] = fileparts(tmpList{i});
    end
    
    if ~isempty(pth)
        %%% sep into sets
        [UQ ia ic] = unique(pth);
        for u = 1:size(UQ,2)
            FileList{end+1} = tmpList(ic==u);
        end
    else
        FileList = {};
    end
    
    
    
end

%{
    FilePath = '/mnt/spaldingimages/';
    FileList = {};
    FileExt = {'tiff','TIF'};
    verbose = 1;
    SET = sdig(FilePath,FileList,FileExt,verbose);
%}