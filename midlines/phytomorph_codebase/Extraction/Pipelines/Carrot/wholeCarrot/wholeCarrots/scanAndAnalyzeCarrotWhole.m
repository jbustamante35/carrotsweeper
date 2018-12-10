function [] = scanAndAnalyzeCarrotWhole(user,auth)

    % whole analysis type
    tissueType = 'wholes';
    % plantType
    plantType = 'carrot';
    % get file list
    [FileList] = ScanAndIssueNewFilesOniRods(user,tissueType,plantType,{'tif','TIF','tiff','nef','NEF'},0);
    
    
    
    
    

    
    %{
    
    % project name
    projectName = 'wholeCarrot';
    % whole analysis type
    tissueType = 'wholes';
    % plantType
    plantType = 'carrot';
    % file ext
    FileExt = {'NEF'};
    % get file list
    [FileList] = ScanAndIssueNewFilesOniRods(user,tissueType,plantType,FileExt,0);
    % transform into URL list
    [FileList] = xform2URL(FileList);
    % geneate the dag
    tmpFileLocation = ['/mnt/scratch1/phytomorph_dev/Deploy/' projectName '/tmpSubmitFiles/'];
    mkdir(tmpFileLocation);
    % output location
    local_oPath = ['/mnt/spaldingdata/nate/mirror_images/carrotData/' user '/return/wholeData/'];
    
    
    
    %[FileList] = xform2URL(FileList);
   
    %}
    numJobs = numel(FileList);
    % remote output location for irods push
    remoteOutputLocation = ['/iplant/home/' user '/#plantType#/return/#tissueType#/'];
    remoteOutputLocation = strrep(remoteOutputLocation,'#plantType#','carrotData');
    remoteOutputLocation = strrep(remoteOutputLocation,'#tissueType#','wholeData');
    CMD = ['imkdir -p ' remoteOutputLocation];
    system(CMD);
    [remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),10*numJobs,'write');
    
    
    
   
    
    % spinup cFlow
    func = cFlow('singleWholeCarrotAnalyze');
    func.setMCRversion('v840');
    func.setMemory(4000);
    numJobs = numel(FileList);
    %numJobs = 10;
    for e = 1:numJobs
        fprintf(['start generating job:' num2str(e) ':' num2str(numJobs) '\n']);
        func(FileList{e},3350,200,40,'./output/',remoteOutputLocation);
        fprintf(['end generating job:' num2str(e) ':' num2str(numJobs) '\n']);
    end
    func.submitDag(auth,700,700);
    
end

%{  
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    scanAndAnalyzeCarrotWhole('turnersarahd',auth);
%}