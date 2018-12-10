classdef cFlow < handle
    
    properties (Access = public)
        jobList;
        jobFunction;
        outputLocation;
        uniqueTimeRandStamp;
        dateString;
        retryValue = 3;
        maxidle = [];
        maxjobs = [];
        maxpost = [];
        
        prescriptFile = '';
        subMemFunc = '';
        tmpFilesLocation;
        submitNodeLocation;
        
        dirMappingsString;
        
        % post script untar line
        mainline0 = '#!/bin/sh';
        utarLine = 'tar xvf "$#N1#" -C $#N2#';
        rmtarLine = 'rm "$#N1#"';
        
        % file dependency lists
        localD = {};
        squidD = {};
        
        % connection to valueDatabase
        valueDatabase = [];
        
        % n route
        n_route_to_shell = 0;
        
        % MCR version
        MCR_version = 'v717';
        
        % memory use
        memUse = '';
        
        % algo name and ver
        algoName = '';
        algoVer = '';
        
        % is GPU
        isGPU = 0;
    end
    
    properties (Constant)
        defaultCompileDirectory = '/mnt/myCompile/';
        %defaultCompileDirectory = '/mnt/scratch1/junkCompile/';
        %defaultCompileDirectory = '/mnt/spaldingdata/nate/inMemCondor/compiledFunctions/';
        %defaultCompileDirectory = '/mnt/snapper/nate/inMemCondor/compiledFunctions/'
        %defaultCompileDirectory = '/mnt/scratch1/nate/inMemCondor/compiledFunctions/';
    end
    
    properties (Access = private)
        dagline_headNode = 'JOB headnode headnode.nothing DONE';
        
        dagline_submitfile = 'JOB #jobName# #jobSubmitfile#';
        
        dagline_vars_filetransferlist = 'VARS #jobName# FileTransferList = "#FileTransferList#"';
        
        dagline_vars_input_line = 'VARS #jobName# argNumber#N# = "#argValue#"';
        dagline_vars_output_line = 'VARS #jobName# argNumber#N# = "#jobName#"';
        dagline_retry = 'RETRY #jobName# #numberRetries#';
        dagline_prescript = 'SCRIPT PRE #jobName# #preScriptName#';
        datline_tar_output = 'SCRIPT POST #jobName# post.sh #jobName#.tar #outputLocation#';
        datline_tar_output_mod = 'SCRIPT POST #jobName# post.sh ';
        dagline_graph_dec = 'PARENT headnode CHILD ';
    end
    
    methods
        function [obj] = cFlow(func,algoName,algoVer,dirMappings,valueDatabase)
            if nargin >= 1
                % generate unique time-rand stamp
                obj.uniqueTimeRandStamp = strrep([num2str(now) num2str(rand(1,1))],'.','');
                obj.dateString = datestr(now);
                obj.subMemFunc = func;
                % set the function name as a string
                setFunctionName(obj,'cFlow_execute');
                % call compile - has built in compile logic
                uniqueDAG = cFlow.compileFunction(func,obj.uniqueTimeRandStamp);
                % set the dag launch directory
                setTempFilesLocation(obj,uniqueDAG);
                % make output location
                obj.outputLocation = cFlow.generateUniqueOutputLocation(func,obj.uniqueTimeRandStamp);
                CMD = ['mkdir -p ' obj.outputLocation];
                system(CMD);
                % make default mapping available for the memory
                obj.addDirectoryMap([cJob.deployed_ouput_vars_location '>' obj.outputLocation]);
                
                
                
            end
            
            
            if nargin == 3
                obj.n_route_to_shell = 1;
                obj.algoName = algoName;
                obj.algoVer = algoVer;
            end
            
            if nargin == 6
                % if the arg is passed in as a string only and not a cell 
                if ~iscell(dirMappings)
                    dirMappings = {dirMappings};
                end
                % add all the directory mappings to the dag
                for e = 1:numel(dirMappings)
                    obj.addDirectoryMap(dirMappings{e});
                end
            end
            
            
            if nargin == 4
                obj.valueDatabase = valueDatabase;
            end
        end
        
        function [] = addDirectoryMap(obj,dirMapString)
            obj.dirMappingsString{end+1} = dirMapString;
        end
        
        function [] = addJob(obj,job)
            obj.jobList{end+1} = job;
        end
        
        function [] = addLocalD(obj,dFile)
           obj.localD{end+1} = dFile; 
        end
        
        function [] = addSquidD(obj,dFile)
            obj.squidD{end+1} = dFile;
        end
        
        function [] = setMCRversion(obj,ver)
            obj.MCR_version = ver;
        end
        
        function [] = setFunctionName(obj,jobFunction)
            obj.jobFunction = jobFunction;
        end
        
        function [] = setOutputLocation(obj,outputLocation)
            obj.outputLocation = outputLocation;
        end
        
        function [] = setTempFilesLocation(obj,tmpLocation)
            obj.tmpFilesLocation = tmpLocation;
        end
        
        function [] = setSubmitNodeLocation(obj,submitNodeLocation)
            obj.submitNodeLocation = submitNodeLocation;
        end
        
        function [] = setMemory(obj,mem)
            obj.memUse = mem;
        end
        
        function [] = setGPU(obj,numGPU)
            obj.isGPU = numGPU;
        end
        
        % this function will render a local copy of the dag files
        function [] = renderDagFile(obj,oFilePath)
            if nargin == 1
                oFilePath = obj.tmpFilesLocation;
            end
            fileID = fopen([oFilePath obj.generate_dagName] ,'w');
            % setup for headnode
            fprintf(fileID,'%s\n',obj.dagline_headNode);            
            %
            tmp_dagline_graph_dec = obj.dagline_graph_dec;
            % renderjobs
            for jb = 1:numel(obj.jobList)
                % assign job name var
                jobName = ['job' num2str(jb)];
                % build up graph dec
                tmp_dagline_graph_dec = [tmp_dagline_graph_dec jobName ' '];
                % assign submitfile var - homogenous waste for now                
                submitFileName = obj.jobList{jb}.generate_submitName;
                
                % setup job name and submit file
                tmp = strrep(obj.dagline_submitfile,'#jobName#',jobName);
                tmp = strrep(tmp,'#jobSubmitfile#',submitFileName);
                fprintf(fileID,'%s\n',tmp);
                
                % setup file transferlist                
                tmp = strrep(obj.dagline_vars_filetransferlist,'#jobName#',jobName);
                tmp = strrep(tmp,'#FileTransferList#',obj.jobList{jb}.getTransferFileList());
                fprintf(fileID,'%s\n',tmp);
                                
                % setup for arguments in
                nargs = obj.jobList{jb}.jobNargin;
                for arg = 1:nargs
                    argValue = obj.jobList{jb}.getArgument(arg);
                    tmp = strrep(obj.dagline_vars_input_line,'#jobName#',jobName);
                    tmp = strrep(tmp,'#N#',num2str(arg));                    
                    tmp = strrep(tmp,'#argValue#',argValue);
                    fprintf(fileID,'%s\n',tmp);
                end
                
                %if isempty(obj.dirMappingsString)
                    % setup for output tar via job name
                    tmp = strrep(obj.dagline_vars_output_line,'#jobName#',jobName);
                    tmp = strrep(tmp,'#N#',num2str(nargs+1)); 
                    fprintf(fileID,'%s\n',tmp);
                %else
                    for e = 1:numel(obj.dirMappingsString)
                        % setup for output tar via job name
                        tmp = strrep(obj.dagline_vars_input_line,'#jobName#',jobName);
                        tmp = strrep(tmp,'#N#',num2str(nargs+1+(e)));
                        tmp = strrep(tmp,'#argValue#',[jobName '_dirMapping' num2str(e)]);
                        fprintf(fileID,'%s\n',tmp);
                    end
                %end
                
                % setup for retry
                tmp = strrep(obj.dagline_retry,'#jobName#',jobName);
                tmp = strrep(tmp,'#numberRetries#',num2str(obj.retryValue));
                fprintf(fileID,'%s\n',tmp);
                
                % if pre script
                if ~isempty(obj.prescriptFile)
                    tmp = strrep(obj.dagline_prescript,'#jobName#',jobName);
                    tmp = strrep(tmp,'#preScriptName#',obj.prescriptFile);
                    fprintf(fileID,'%s\n',tmp);    
                end
                
                if isempty(obj.dirMappingsString)
                    % setup for untar output in post script
                    tmp = strrep(obj.datline_tar_output,'#jobName#',jobName);
                    tmp = strrep(tmp,'#outputLocation#',obj.outputLocation);
                    fprintf(fileID,'%s\n',tmp);
                else
                    tmpTarName = strrep(obj.datline_tar_output_mod,'#jobName#',jobName);
                    for e = 1:numel(obj.dirMappingsString)
                        fidx = strfind(obj.dirMappingsString{e},'>');
                        tmpTarName = [tmpTarName jobName '_dirMapping' num2str(e) '.tar' ' ' obj.dirMappingsString{e}(fidx(1)+1:end) ' '];
                        
                    end
                    fprintf(fileID,'%s\n',tmpTarName);
                end
            end
            
            if ~isempty(obj.maxidle)
                fprintf(fileID,'%s\n',['maxidle=' num2str(obj.maxidle)]);
            end
             
            if ~isempty(obj.maxjobs)
                fprintf(fileID,'%s\n',['maxjobs=' num2str(obj.maxjobs)]);
            end
             
            if ~isempty(obj.maxpost)
                fprintf(fileID,'%s\n',['maxpost=' num2str(obj.maxpost)]);
            end
            
            fprintf(fileID,'%s\n',tmp_dagline_graph_dec);
            fclose(fileID);
        end
        
        function [dagName] = generate_dagName(obj)
            dagName = [obj.jobFunction '.dag'];
        end
        
        function [scpFileList] = generate_scpFileList(obj)
            scpFileList{1} = [obj.tmpFilesLocation obj.generate_dagName()];
            scpFileList{2} = [obj.tmpFilesLocation obj.jobList{1}.generate_submitName()];
            scpFileList{3} = [obj.tmpFilesLocation obj.jobList{1}.generate_exeName()];
            scpFileList{4} = [obj.tmpFilesLocation obj.jobFunction];
            scpFileList{5} = [obj.tmpFilesLocation 'run_' obj.jobFunction '.sh'];
            CMD = ['sed -i ''s/LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}\/runtime\/glnxa64:${MCRROOT}\/sys\/opengl\/lib\/glnxa64;/LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}\/sys\/opengl\/lib\/glnxa64:${PWD}\/SS:${PWD}\/lcms\/lib;/g'' ' scpFileList{5}];
            system(CMD);
            scpFileList{6} = [obj.defaultCompileDirectory 'clear.sh'];
            scpFileList{7} = [obj.tmpFilesLocation 'post.sh'];
            if ~isempty(obj.prescriptFile)
                scpFileList{6} = [obj.tmpFilesLocation obj.prescriptFile];
            end
        end
        
        function [] = submitDag(obj,icommands_auth,varargin)
            
            if nargin >= 3
                maxidle = varargin{1};
            end
            if nargin >= 4
                maxpost = varargin{2};
            end
            remote_DAG_location = [obj.jobFunction filesep obj.uniqueTimeRandStamp];
            
            % 
            if isempty(obj.dirMappingsString)
                obj.jobList{1}.generate_submitFilesForDag(icommands_auth);
            else
                obj.jobList{1}.generate_submitFilesForDag(icommands_auth,obj.dirMappingsString);
            end
            
            
            obj.renderDagFile();
            obj.generatePostScript();
            scpList = obj.generate_scpFileList();
            
            
            
            
            
            dirCMD_logs_out = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stdout/'''];
            dirCMD_logs_err = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stderr/'''];
            dirCMD_output = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/output/'''];
            [status result] = system(strrep(dirCMD_logs_out,'#directory#',remote_DAG_location));
            [status result] = system(strrep(dirCMD_logs_err,'#directory#',remote_DAG_location));
            [status result] = system(strrep(dirCMD_output,'#directory#',remote_DAG_location));
            dirCMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir /home/nate/condorFunctions/#directory#/'''];
            [status result] = system(strrep(dirCMD,'#directory#',remote_DAG_location));
            CMD = 'scp -P 50118 #srcfile# nate@128.104.98.118:/home/nate/condorFunctions/#directory#/#desfile#';
            CMD = strrep(CMD,'#directory#',remote_DAG_location);
            for f = 1:numel(scpList)
                [pth nm ext] = fileparts(scpList{f});
                tCMD = strrep(CMD,'#desfile#',[nm ext]);
                tCMD = strrep(tCMD,'#srcfile#',scpList{f});
                [status result] = system(tCMD);
            end
           
            
            
            
            % submit the job dag
            dagName = obj.generate_dagName();
            CMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'cd /home/nate/condorFunctions/#directory#/; condor_submit_dag -maxidle ' num2str(maxidle) ' -maxpost ' num2str(maxpost) ' ' dagName ''''];
            CMD = strrep(CMD,'#directory#',remote_DAG_location);
            system(CMD,'-echo');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % condor launch - END
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function [] = generateSleepScript(obj,sleepTime,units)
            sleepFile = [obj.tmpFilesLocation 'sleepScript.sh'];
            fileID = fopen(sleepFile,'w');
            fprintf(fileID,'%s\n',['sleep ' num2str(sleepTime) units]);
            fclose(fileID);
            obj.prescriptFile = 'sleepScript.sh';
        end
        
        function [varargout] = subsref(obj,s)
            if strcmp(s(1).type,'()')
                
                
                
                
                tmpJob = cJob();
                % set GPU
                tmpJob.setAsMemoryJob(obj.uniqueTimeRandStamp,obj.subMemFunc);
                tmpJob.setNumberofArgs(numel(s.subs)-obj.n_route_to_shell);
                tmpJob.setNumberofOArgs(nargout);
                
                
                %[obj] = fargList(s.subs,obj.valueDatabase,tmpJob.fullMatLocation);
                
                
                %tmpJob.fullMatLocation
                for e = 1:(numel(s.subs)-obj.n_route_to_shell)
                    argInputString{e} = tmpJob.setArg(s.subs{e+obj.n_route_to_shell},e);
                end
                for e = 1:nargout
                    % add the output directory to the string for var name
                    %varargout{e} = [obj.outputLocation 'output' filesep tmpJob.matFileName '@out' num2str(e)];
                    varargout{e} = [obj.outputLocation tmpJob.matFileName '@out' num2str(e)];
                end
                save(tmpJob.fullMatLocation,'tmpJob','-v7.3','-append');
                
                
                
                
                % make a cJob which will wrap the cJob - why?
                %wrapJob = cJob();
                wrapJob = cJob(obj.algoName,obj.algoVer);
                wrapJob.isGPU = obj.isGPU;
                if ~isempty(obj.memUse)
                    wrapJob.setMemoryUse(obj.memUse);
                end
                wrapJob.changeMCRfile(obj.MCR_version);
                for e = 1:numel(obj.localD)
                    wrapJob.addFile(obj.localD{e});
                end
                for e = 1:numel(obj.squidD)
                    wrapJob.addSquidFile(obj.squidD{e});
                end
                wrapJob.setTempFilesLocation(obj.tmpFilesLocation);
                wrapJob.addFile(tmpJob.fullMatLocation);
                wrapJob.setFunctionName('cFlow_execute');    
                wrapJob.setNumberofArgs(1+obj.n_route_to_shell);
                wrapJob.setArgument(tmpJob.matFileName,1);
                % direct value vs pointer
                if obj.n_route_to_shell ~= 0
                    
                    
                    pointerValue = s.subs{1};
                    wrapJob.setArgument(pointerValue,2);
                    
                    
                    % added for miron
                    %wrapJob.addFile(s.subs{1});
                    
                end
                %{
                if isempty(obj.dirMappingsString)
                    wrapJob.generate_submitFilesForDag(icommands_auth);
                else
                    wrapJob.generate_submitFilesForDag(icommands_auth,obj.dirMappingsString);
                end
                %}
                obj.addJob(wrapJob);
            else
                [varargout{1:nargout}] = builtin('subsref',obj,s);
            end
        end
        
        function [] = generatePostScript(obj)
            remote_DAG_location = [obj.jobFunction filesep obj.uniqueTimeRandStamp];
            fileID = fopen([obj.tmpFilesLocation 'post.sh'],'w');
            % setup for post
            fprintf(fileID,'%s\n',obj.mainline0);
            for e = 1:numel(obj.dirMappingsString)
                fidx = strfind(obj.dirMappingsString{e},'>');
                s1 = num2str((e-1)*2+1);
                s2 = num2str((e-1)*2+2);
                utarLine = strrep(obj.utarLine,'#N1#',[s1]);
                utarLine = strrep(utarLine,'#N2#',[s2  ' ' obj.dirMappingsString{e}(1:(fidx(1)-1)) '/* --strip-components=1']);
                % setup for untarline
                fprintf(fileID,'%s\n',utarLine);
                utarLine = strrep(utarLine,'#N2#',s2);
                % setup for rmtarline
                rmtarLine = strrep(obj.rmtarLine,'#N1#',s1);
                fprintf(fileID,'%s\n',rmtarLine);
            end
        end
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compile 
        % inputs: func: string that will need to be compiled
        %         uniqueTimeRandStamp: the dag for this evaluation of the function  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [uniqueEvalDirectory] = compileFunction(func,uniqueTimeRandStamp)
            tmpCompileDirectory = '/mnt/scratch1/junkCompile/';
            
            % get the location of the function that is being used
            pathToFunction = which(func);
            % get the default compile directory
            compile_directory = cFlow.defaultCompileDirectory;
            % make the function compile directory
            functionDirectory = [compile_directory func filesep];
            mkdir(functionDirectory);
            % make the uniqueTimeRandStamp directory
            uniqueEvalDirectory = [functionDirectory 'DAG' filesep uniqueTimeRandStamp filesep];
            CMD = ['mkdir -p ' uniqueEvalDirectory];
            system(CMD);
            % make the backup directory
            functionBackupDirectory = [functionDirectory 'backUp' filesep];
            mkdir(functionBackupDirectory);
            % make the function file
            functionFile = [functionDirectory func '.m'];
            % set the default compile flat to false
            compileFlag = 0;
            % check to see if the functionFile already is present
            if exist(functionFile)
                % if the current version in the cFlow repo is different
                % from the known version, then backup the current, and
                % compile
                CMD = ['diff ' functionFile ' ' pathToFunction];
                [status, result] = system(CMD);
                if ~isempty(result)
                    % make backup of function
                    dateStamp = datestr(now);
                    dateStamp = strrep(dateStamp,':','_');
                    dateStamp = strrep(dateStamp,'-','_');
                    dateStamp = strrep(dateStamp,' ','_');
                    [p,n,e] = fileparts(functionFile);
                    CMD = ['cp ' functionFile ' ' functionBackupDirectory n '_' dateStamp '.m'];
                    system(CMD);
                    % copy in new function
                    CMD = ['cp ' pathToFunction ' ' functionFile];
                    system(CMD);
                    compileFlag = 1;
                end
            else
                CMD = ['cp ' pathToFunction ' ' functionFile];
                system(CMD);
                compileFlag = 1;
            end
            if compileFlag
                tmpCompileDirectory = [tmpCompileDirectory filesep uniqueTimeRandStamp filesep];
                CMD = ['mkdir -p ' tmpCompileDirectory];
                system(CMD)
                %CMD = ['mcc -d ' functionDirectory ' -m -v -R -singleCompThread ' functionFile];
                %CMD = ['mcc -d ' tmpCompileDirectory ' -a ' functionFile ' -a  gmdistribution.m -a cJob.m -a im2single.m -m -v -R -singleCompThread -R -nojvm cFlow_execute.m'];
                %CMD = ['mcc -d ' tmpCompileDirectory ' -m -v -R -singleCompThread -a ' functionFile ' -a SeriesNetwork.m -a extractSingleBugEye.m -a condorDeploy_ver0.m -a partialFunction.m -a  gmdistribution.m -a cJob.m -a im2single.m -m -v -R -singleCompThread cFlow_execute.m'];
                CMD = ['mcc -d ' tmpCompileDirectory ' -m -v -R -singleCompThread -a ' functionFile ' -a applyAllLayers.m -a SeriesNetwork.m -a ClassificationDiscriminant.m -a partialFunction.m -a current_sorghum_network.m -a generalizeLoader.m -a generalizeFeatureExtractor.m -a gmdistribution.m -a cJob.m -a im2single.m -m -v -R -singleCompThread cFlow_execute.m'];
                
                eval(CMD);
                
                % copy the compiled function and its nessary scripts into the
                % evalution directory
                filesToCopy = dir(tmpCompileDirectory);
                filesToCopy([filesToCopy.isdir]) = [];
                for e = 1:numel(filesToCopy)
                    sourceFile = [tmpCompileDirectory filesToCopy(e).name];
                    targetFile = [functionDirectory filesToCopy(e).name];
                    CMD = ['cp ' sourceFile ' ' targetFile];
                    system(CMD);
                end
            end
            
            %{
            uniqueEvalDirectory = strrep(uniqueEvalDirectory,'/mnt/scratch1/junkCompile/','/mnt/myCompile/');
            CMD = ['mkdir -p ' uniqueEvalDirectory];
            system(CMD);
            %}
            
            % copy the compiled function and its nessary scripts into the
            % evalution directory
            filesToCopy = dir(functionDirectory);
            filesToCopy([filesToCopy.isdir]) = [];
            for e = 1:numel(filesToCopy)
                sourceFile = [functionDirectory filesToCopy(e).name];
                targetFile = [uniqueEvalDirectory filesToCopy(e).name];
                CMD = ['cp ' sourceFile ' ' targetFile];
                system(CMD);
            end
            
            
        end
        
        function [uniqueEvalDirectory] = generateUniqueCompileLocation(func,uniqueTimeRandStamp)
            % get the default compile directory
            compile_directory = cFlow.defaultCompileDirectory;
            % make the function compile directory
            functionDirectory = [compile_directory func filesep];
            % make the uniqueTimeRandStamp directory
            uniqueEvalDirectory = [functionDirectory 'DAG' filesep uniqueTimeRandStamp filesep];
        end
        
        function [uniqueOutputLocation] = generateUniqueOutputLocation(func,uniqueTimeRandStamp)
            [uniqueEvalDirectory] = cFlow.generateUniqueCompileLocation(func,uniqueTimeRandStamp);
            uniqueOutputLocation = [uniqueEvalDirectory 'functionOutputs' filesep];
        end
        
        function [] = persistFunctionCall()
            
        end
    end
end

%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple test - many-inputs > one-output
    func = cFlow('testCondorFunction');
    res = func(1,2);
    func.submitDag(50,50);
    o = cFlowLoader(res);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple test - many-inputs > many-outputs
    func = cFlow('testCondorFunction');
    [res1,res2] = func(1,2);
    func.submitDag(50,50);
    [o1,o2] = cFlowLoader(res1,res2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple test - loop over inputs
    func = cFlow('testCondorFunction');
    for e = 1:10
        [resM1{e} resM2{e}] = func(e,e+1);
    end
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    func.submitDag(auth,50,50);
    for e = 1:10
        [o1M{e} o2M{e}] = cFlowLoader(resM1{e},resM2{e});
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test directory mapping(s)
    func = cFlow('testCondorFunction',{'saveTo>/mnt/spaldingdata/nate/saveTo/','saveTo2>/mnt/spaldingdata/nate/saveTo2/'});
    res = func(1,2,'./saveTo/','./saveTo2/');
    func.submitDag(50,50);
    o = cFlowLoader(res);

    




   
%}