function [FileList rawFileList] = ScanAndIssueNewFilesOniRods(user,tissueType,plantType,fileExt,removedProcessed)
    
    if nargin == 2
        plantType = 'maize';
        fileExt = {'tif','TIF','tiff','nef'};
    end
    
    if nargin == 2 & strcmp(tissueType,'seedlings')
        fileExt = {'nef'};
    end

    if nargin <= 4
        removedProcessed = 1;
    end
    
    % scan a fuse mounted directory for maize data
    [FileList] = scanForImagesOnIrods(user,plantType,tissueType,fileExt);
    %FileList = FileList(1:50);
    
    %CMD = uncLog(FileList,'add',[plantType '-' tissueType],'1.0',[],[],1);
    %{
    % remove the processed files
    if removedProcessed
        [FileList] = removeProcessedFiles(FileList,user,tissueType,plantType);
    end
    
    
    % change the local file name to irods file name and issue tickets
    baseString = ['/iplant/home/' user '/$plantTypeData/'];
    baseString = strrep(baseString,'$plantType',plantType);
    [FileList] = fuse2irods2(FileList,'/home/nate/iplant/',baseString);
    %}
    % issue bulk shared tickets
    fprintf(['Starting ticket issue \n']);
    rawFileList = FileList;
    [FileList] = issueBulkTicket(FileList);
    fprintf(['Ending ticket issue \n']);
    
end