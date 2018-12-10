function [] = scanAndAnalyzeAutoKinematics(user)
    % project name
    projectName = 'straightKinematics';
    % whole analysis type
    tissueType = 'rootStraight';
    % plantType
    plantType = 'arabidopsis';
    % file ext
    FileExt = {'tif','TIF','png'};
    % get file list
    [FileList] = ScanAndIssueNewFilesOniRods(user,tissueType,plantType,FileExt);
    % transform into URL list
    %[FileList] = xform2URL(FileList);
    % order into movie sets
    [FileList] = orderFrom_gdig(FileList,{});
    % geneate the dag
    tmpFileLocation = ['/mnt/scratch1/phytomorph_dev/Deploy/' projectName '/tmpSubmitFiles/'];
    CMD = ['mkdir -p ' tmpFileLocation];
    system(CMD);
    % output location
    local_oPath = ['/mnt/spaldingdata/nate/mirror_images/arabidopsisData/' user '/return/straightKinematicsQ1/'];
    CMD = ['mkdir -p ' local_oPath];
    system(CMD);
    % spinup cFlow
    func = cFlow('singleKinematicsStack',{['output>' local_oPath]});
    %func.setMCRversion('v840');
    for e = 1:numel(FileList)
        func(FileList{e},0,'./output/',1,10);
    end
    func.submitDag(150,150);
end

%{
    scanAndAnalyzeAutoKinematics('monshausenlab');
    scanAndAnalyzeAutoKinematics('trieupham10');
    oPath = ['/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/return/straightKinematicsW'];
    singleKinematicsStack(FileList{1},0,oPath,1,10);


    
FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/straight/hypoosmotic/';
FileList = {};
FileExt = {'tif'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
%}