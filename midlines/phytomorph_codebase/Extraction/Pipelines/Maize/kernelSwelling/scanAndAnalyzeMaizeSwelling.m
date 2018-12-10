function [] = scanAndAnalyzeMaizeSwelling(user,local)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scan for new images to analyze - START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~local
        % get file list, transform file list for URL, issue tickets for write
        [FileList] = ScanAndIssueNewFilesOniRods(user,'kernelSwellData','maize',{'tif','TIF','tiff','nef'},0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% scan for new images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        FileList = {};
        FileExt = {'tiff','TIF','tif'};
        verbose = 1;
        SET = sdig(inFilePath,FileList,FileExt,verbose);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% sort SET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        N = [];
        for img = 1:numel(SET{e})
            [p n ex] = fileparts(SET{e}{img});
            N(img) = str2num(n);
        end
        [N sidx] = sort(N);
        SET{e} = SET{e}(sidx);
    end
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% remove those which are in the junk folder
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        fidx = strfind(SET{e}{1},'junk');
        if ~isempty(fidx)
            rmidx(e) = 1;
        else
            rmidx(e) = 0;
        end
    end
    SET(find(rmidx)) = [];    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% find sets with more than 250
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        numImages(e) = numel(SET{e});
    end
    rmidx = comFunc(numImages,expectedImageNumber);
    SET(rmidx) = [];
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% generate output file names
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        [swell{e} area{e} para{e} err{e} fit{e}] = generateOutFileBase(SET{e}{1},oPath,1);
    end
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% remove sets which have data in the oPath
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% remove compiled results
    oFileList = {};
    FileExt = {'csv'};
    verbose = 1;
    oSET = gdig(oPath,FileList,FileExt,verbose);
    rmidx = [];
    for e = 1:numel(oSET)
        fidx = strfind(oSET{e},'compiled_results');
        if ~isempty(fidx)
            rmidx(e) = 1;
        else
            rmidx(e) = 0;
        end
    end
    oSET(find(rmidx)) = [];
    %}
    %{
    %%% look for data which already is run through algo
    rmidx = zeros(numel(SET),1);    
    for e = 1:numel(oSET)
        fidx1 = strfind(oSET{e},filesep);
        fidx2 = strfind(oSET{e},'--');
        snipFile = oSET{e}(fidx1(end)+1:fidx2(end)-1);
        for i = 1:numel(swell)
            fidx1 = strfind(swell{i},filesep);
            fidx2 = strfind(swell{i},'--');
            isnipFile = swell{i}(fidx1(end)+1:fidx2(end)-1);
            if strcmp(snipFile,isnipFile)
                rmidx(i) = 1;
            end
        end
    end
    
    SET(find(rmidx)) = [];
    swell(find(rmidx)) = [];
    area(find(rmidx)) = [];
    para(find(rmidx)) = [];
    err(find(rmidx)) = [];
    fit(find(rmidx)) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scan for new images to analyze - END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    %{
    
    % generate dag
    tmpFileLocation = '/mnt/scratch1/phytoM/flashProjects/maizePipeline/maizeKernel/tmpSubmitFiles/';
    dag = epfod();
    dag.setFunctionName('singleKernelImage');
    dag.setOutputLocation(['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/kernelData/']);
    dag.setTempFilesLocation(tmpFileLocation);
    % add jobs to dag for each image - create and add job to dag
    numJobs = numel(FileList);
    for e = 1:numJobs
        [pth,nm,ex] = fileparts(FileList{e});
        % create job
        job = cJob();
        job.setTempFilesLocation(tmpFileLocation);
        job.setFunctionName('singleKernelImage');    
        job.setNumberofArgs(4);
        job.setArgument(FileList{e},1);        
        job.setArgument('./output/',2); 
        job.setArgument('1',3);
        job.setArgument('1',4); 
        % add job to dag
        dag.addJob(job);
        job.generate_submitFilesForDag();
    end
    % submit dag
    dag.submitDag(50,50);
    %{
        singleKernelImage(FileList{73},'./output/',1,1);
    %}
    %}
end

%{   
    
    scanAndAnalyzeMaizeSwelling('jgustin');
%}