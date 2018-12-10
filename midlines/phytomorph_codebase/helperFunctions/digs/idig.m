function [FileList] = idig(FilePath,FileExt)    
    cmd = ['iquest --no-page "SELECT DATA_ID, DATA_NAME, COLL_NAME WHERE COLL_NAME like ''' FilePath '%''"'];
    [status query] = system(cmd);
    if ~strcmp(FilePath,'/')
        FilePath = [FilePath '/'];
    end
    verbose = 1;
    FileList = [];
    cdir = mydir(FilePath);
    [pth,nm,ext] = fileparts(FilePath);
    FileList = sdir(cdir,FileExt,pth);
end