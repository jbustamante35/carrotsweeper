function [] = scanForImages(in)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % look for image data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %%%% read the csv file containing the paths to scan
    toScan = readtext(in.inScanFile);    
    
    %%%% setup file scanner - change to scan for iplant later    
    FileExt = {'tiff','TIF','tif','JPG'};           % the file extensions to scan for    
    verbose = 1;                                    % if verbose
    
    fprintf(['Launch arabidopsis pipeline:\n']);
    for e = 1:size(toScan)
        fprintf(['[' num2str(e) ']-' toScan{e} '\n']);
    end
    evalResponse = input('Please select pipeline for analysis:');
    
    %%%% scan each path for image sets
    masterSet = {};
    %for s = 1:size(toScan,1)        
    for s = evalResponse;
        % perform scan
        fileSet = sdig(toScan{s},{},FileExt,verbose);
        % sort the image stacks
        fileSet = sortSets(fileSet);
        masterSet = [masterSet,fileSet];
    end
    
    outputLocation = ['/mnt/snapper/pipeLines/arabidopsis/return/output/'];
    %%%% setup file scanner - change to scan for iplant later    
    FileExt = {'mat'};         % the file extensions to scan for    
    verbose = 1;                            % if verbose
    resultsFileList = gdig(outputLocation,{},FileExt,verbose);
    toRM = [];
    for s = 1:numel(masterSet)
        [tp tn te] = fileparts(masterSet{s}{1});
        tmp = strrep(tp,filesep,'--');
        tmp = [outputLocation tmp '.mat'];
        if sum(strcmp(tmp,resultsFileList))
           toRM = [toRM s]; 
        end
    end
    masterSet(toRM) = [];
    
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condor launch - START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate the dag 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate dag
    tmpFileLocation = '/mnt/scratch1/phytoM/flashProjects/arabidopsis_auto/';
    dag = epfod();
    dag.setFunctionName('isolateRoots_overStack');
    %dateString = strrep(strrep(strrep(datestr(clock, 'yyyy-mm-dd HH:MM:SS'),' ','_'),':','_'),'-','_');
    outputLocation = ['/mnt/snapper/share/pipeLines/arabidopsis/return/'];
    mkdir(outputLocation);
    dag.setOutputLocation(outputLocation);
    dag.setTempFilesLocation(tmpFileLocation);

    numJobs = numel(masterSet);
    %numJobs = 100;
    for e = 1:numJobs

        %%%%%%%%%%%%%%%%%%%%%
        % stack - the file names of the images
        % outPath - the location to write the results to
        % toWrite - flag for writing csv files    
        % NP - number of points for tip angle measurement
        % NPK - number of points to measure kurvature over
        % SNIP - number of points to SNIP from midline
        % set job arguments
        toWrite = '1';
        NP = '20';
        NPK = '20';
        SNIP = '20';
        disp = '0';
        
        %%%%%%%%%%%%%%%%%%%%%
        % file(s) into parts
        clear condorLocalNames
        for f = 1:numel(masterSet{e})        
            [pth,condorLocalNames{f},ex] = fileparts(masterSet{e}{f});
            condorLocalNames{f} = [condorLocalNames{f} ex];
        end
        condorLocalNames = cellStringConvert(condorLocalNames);

        %%%%%%%%%%%%%%%%%%%%%
        % generate output name
        pth = strrep(pth,filesep,'SLASH');
        saveName = pth;
        saveName = strrep(saveName,' ','SPACE');

        %%%%%%%%%%%%%%%%%%%%%
        % create job
        job = cJob();
        job.setOSG(1);
        job.setFlock(1);
        job.setTempFilesLocation(tmpFileLocation);
        job.setFunctionName('isolateRoots_overStack');    
        job.setNumberofArgs(9);

        %%%%%%%%%%%%%%%%%%%%%
        % add image file(s) to job
        for f = 1:numel(masterSet{e})
            job.addFile(masterSet{e}{f});
        end

        %%%%%%%%%%%%%%%%%%%%%
        % set arguments
        job.setArgument(condorLocalNames,1);        
        job.setArgument('./output/',2);
        job.setArgument(saveName,3);
        job.setArgument(toWrite,4);
        job.setArgument(NP,5);
        job.setArgument(NPK,6);
        job.setArgument(SNIP,7);
        job.setArgument(disp,8);
        job.setArgument('-1',9);

        % add job to dag
        dag.addJob(job);
        job.generate_submitFilesForDag();
    end

    if numJobs > 0

        dag.renderDagFile();
        scpList = dag.generate_scpFileList();
        dirCMD_logs_out = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stdout/'''];
        dirCMD_logs_err = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stderr/'''];
        dirCMD_output = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/output/'''];
        [status result] = system(strrep(dirCMD_logs_out,'#directory#',dag.jobFunction));
        [status result] = system(strrep(dirCMD_logs_err,'#directory#',dag.jobFunction));
        [status result] = system(strrep(dirCMD_output,'#directory#',dag.jobFunction));
        dirCMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir /home/nate/condorFunctions/#directory#/'''];
        [status result] = system(strrep(dirCMD,'#directory#',dag.jobFunction));
        CMD = 'scp -P 50118 #srcfile# nate@128.104.98.118:/home/nate/condorFunctions/#directory#/#desfile#';
        CMD = strrep(CMD,'#directory#',dag.jobFunction);
        for f = 1:numel(scpList)
            [pth nm ext] = fileparts(scpList{f});
            tCMD = strrep(CMD,'#desfile#',[nm ext]);
            tCMD = strrep(tCMD,'#srcfile#',scpList{f});
            [status result] = system(tCMD);
        end

        % submit the job dag
        dagName = dag.generate_dagName();
        CMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'cd /home/nate/condorFunctions/#directory#/; condor_submit_dag -maxidle 75 -maxjobs 75 -maxpost 50 ' dagName ''''];
        CMD = strrep(CMD,'#directory#',dag.jobFunction);
        system(CMD,'-echo');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % condor launch - END
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    
    %{
    % run local
    % input list for isolateRoots_overStack
    % (stack,outPath,outName,toWrite,NP,NPK,SNIP,disp,numToProcess)
    % stack - the file names of the images    
    % outPath - the location to write the results to
    % outName - the name of the files to save the results as
    % toWrite - flag for writing csv files    
    % NP - number of points for tip angle measurement
    % NPK - number of points to measure kurvature over
    % SNIP - number of points to SNIP from midline
    % disp - to display
    % numImagesToProcess - if less than zero then process all
    %%%% run analysis on each set
    for s = 1:numel(masterSet)
    try
        tm = clock;
        
        out{s} = isolateRoots_overStack(fileSet{s},'./outputSSu1/',num2str(s),1,20,20,20,1,-1);
        etm = etime(clock,tm);
    catch
    end
    end
    %}
end



function [fileSet] = sortSets(fileSet)    
    rmidx = [];
    for e = 1:numel(fileSet)
        try
            NAME = [];
            for f = 1:numel(fileSet{e})
                [pth,nm,ext] = fileparts(fileSet{e}{f});
                NAME(f) = str2num(nm);
            end
            [J sidx] = sort(NAME);
            fileSet{e} = fileSet{e}(sidx);            
        catch ME
            fprintf('error sorting set \n');
            fprintf(['removing set' pth]);
            rmidx = [rmidx e];
        end
    end
    fileSet(rmidx) = [];
end



%{
    in.inScanFile = '/mnt/snapper/pipeLines/arabidopsis/scanPaths.csv';
    in.processedFile = '/mnt/snapper/pipeLines/arabidopsis/processed.csv';
    scanForImages(in);
    % note recompile fix for DC.m where variance in surface curvature is
    % measured rather than variance in image intensity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compile_directory = '/mnt/scratch1/phytoM/flashProjects/arabidopsis_auto/';
    CMD = ['mcc -d ' compile_directory ' -m -v -R -singleCompThread isolateRoots_overStack.m'];
    eval(CMD);

%}