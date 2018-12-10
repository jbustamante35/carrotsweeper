function [FileList] = pdig(FilePath,FileExt)
    FilePath = attachSlash(FilePath);
    verbose = 1;
    FileList = [];
    cdir = mydir(FilePath);
    [pth,nm,ext] = fileparts(FilePath);
    %FileList = sdir(cdir,FileExt,'');
    FileList = struct2list(cdir,'name');
end