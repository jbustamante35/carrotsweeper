function [] = HTmecka(user,algorithm)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                HTmecka.m is create jobs for condor to run Maize Ear Cob Kernel Analysis. A user is 
                to choose one of three algorithm to use and provide user for condor.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                mecka.m, ScanAndIssueNewFilesOniRods.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                user:           The user name for condor. (?)
                algorithm:      The argorithm to use. 'c' for singleCobImage.m, 'e' for singleEarImage.m, and 'k' for singleKernelImage.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    tmpFileLocation = '/mnt/scratch1/maizePipeline/mecka/tmpSubmitFiles/';
    switch algorithm
        case 'c'
            analysisType = 'cobs';
            memREQ = '4000';
            algorithmFlag = 'c';
            numberOfObjects = '3';
            imageRES = '1200';
            localOutputLocation = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/cobData/'];
        case 'e'
            analysisType = 'ears';
            memREQ = '4000';
            algorithmFlag = 'e';
            numberOfObjects = '3';
            imageRES = '1200';
            localOutputLocation = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/earData/'];
        case 'k'
            analysisType = 'kernels';
            memREQ = '4000';
            algorithmFlag = 'k';
            numberOfObjects = [];
            imageRES = '1200';
            localOutputLocation = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/kernelData/'];
    end
    
    
    
    
    % get file list
    [FileList] = ScanAndIssueNewFilesOniRods(user,analysisType);
    % geneate the dag
    dag = epfod();
    dag.setFunctionName('mecka');
    dag.setOutputLocation(localOutputLocation);
    dag.setTempFilesLocation(tmpFileLocation);
    numJobs = numel(FileList);
    % add jobs to dag for each image - create and add job to dag
    for e = 1:numJobs
        [pth,nm,ex] = fileparts(FileList{e});
        % create job
        job = cJob();
        job.requirements.memory = {'=' memREQ};
        job.setTempFilesLocation(tmpFileLocation);
        job.setFunctionName('mecka');
        job.setNumberofArgs(8);
        job.setArgument(algorithmFlag,1);
        job.setArgument([FileList{e}],2);        
        job.setArgument(numberOfObjects,3);
        job.setArgument('./output/',4);
        job.setArgument('1',5);
        job.setArgument('1',6);
        job.setArgument(imageRES,7);
        job.setArgument('1',8)
        
        % add job to dag
        dag.addJob(job);
        job.generate_submitFilesForDag();
    end
    % submit dag
    dag.submitDag(50,50);
end

%{
    HTmecka('gxe','c');
    HTmecka('gxe','e');
    HTmecka('gxe','k');
%}