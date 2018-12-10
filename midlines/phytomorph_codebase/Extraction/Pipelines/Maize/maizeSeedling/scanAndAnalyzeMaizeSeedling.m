function [] = scanAndAnalyzeMaizeSeedling(user,auth)
    analysisType = 'seedlings';
    % get file list
    [FileList] = ScanAndIssueNewFilesOniRods(user,analysisType,'maize',{'tif','TIF','tiff','nef'},0);
    
    remoteOutputLocation = ['/iplant/home/' user '/#plantType#/return/#tissueType#/'];
    remoteOutputLocation = strrep(remoteOutputLocation,'#plantType#','maizeData');
    remoteOutputLocation = strrep(remoteOutputLocation,'#tissueType#','seedlingData');
    
    numJobs = numel(FileList);
    [remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numJobs,'write');
    
    
    %{
    % geneate the dag
    tmpFileLocation = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/tmpSubmitFiles/';
    dag = epfod();    
    dag.setFunctionName('singleSeedlingImage');
    dag.setOutputLocation(['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/seedlingData/']);
    dag.setTempFilesLocation(tmpFileLocation);
    %}
    
    
    %numJobs = 10;
    func = cFlow('singleSeedlingImage');
    func.setMCRversion('v840');
    func.setMemory(2000);
    for e = 1:numJobs
        fprintf(['start generating job:' num2str(e) ':' num2str(numJobs) '\n']);
        func(FileList{e},100,5,400,100,4,20,'./output/',remoteOutputLocation);
        fprintf(['end generating job:' num2str(e) ':' num2str(numJobs) '\n']);
    end
    func.submitDag(auth,500,500);
    
    
    %{
    %numJobs = 10;
    %numJobs = 100;
    % add jobs to dag for each image - create and add job to dag
    for e = 1:numJobs        
        % create job
        job = cJob();
        job.addFile('/mnt/spaldingdata/nate/dcraw');
        job.addSquidFile('core-3.2.1.jar');
        job.addSquidFile('javase-3.2.1.jar');
        job.changeMCRfile('v840');
        job.requirements.memory = {'=' '2000'};
        job.setTempFilesLocation(tmpFileLocation);
        job.setFunctionName('singleSeedlingImage');    
        job.setNumberofArgs(7);        
        job.setArgument([FileList{e}],1);
        job.setArgument('100',2);
        job.setArgument('5',3);
        job.setArgument('120',4);
        job.setArgument('100',5);
        job.setArgument('4',6);
        job.setArgument('./output/',7);
        % add job to dag
        dag.addJob(job);
        job.generate_submitFilesForDag();
    end
    if numJobs ~= 0
        % submit dag
        dag.submitDag(150,150);
    end
    %}
end

%{  
    scanAndAnalyzeMaizeSeedling('hirsc213',auth);

    oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/return/seedlingData/output/';
    parfor e = 1:numel(FileList)
        singleSeedlingImage(FileList{e},100,5,100,100,4,oPath);
    end
   
    singleSeedlingImage(FileList{end},50,5,100,100,4,'/mnt/scratch1/phytoM/output/junk/');
%}