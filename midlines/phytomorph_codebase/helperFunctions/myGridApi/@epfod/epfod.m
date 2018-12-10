classdef epfod < handle
    
    properties (Access = public)
        jobList;
        jobFunction;
        outputLocation;
        
        retryValue = 3;
        maxidle = [];
        maxjobs = [];
        maxpost = [];
        
        prescriptFile = '';
        
        tmpFilesLocation;
        submitNodeLocation;
    end
    
    properties (Access = private)
        dagline_headNode = 'JOB headnode headnode.nothing DONE';
        
        dagline_submitfile = 'JOB #jobName# #jobSubmitfile#';
        
        dagline_vars_filetransferlist = 'VARS #jobName# FileTransferList = "#FileTransferList#"';
        
        dagline_vars_input_line = 'VARS #jobName# argNumber#N# = "#argValue#"';
        dagline_vars_output_line = 'VARS #jobName# argNumber#N# = "#jobName#"';
            
        dagline_retry = 'RETRY #jobName# #numberRetries#';
        dagline_prescript = 'SCRIPT PRE #jobName# #preScriptName#';
        datline_tar_output = 'SCRIPT POST #jobName# post.sh "#jobName#".tar #outputLocation#';
        
        dagline_graph_dec = 'PARENT headnode CHILD ';
    end
    
    methods
        function [obj] = epfod()
            
        end
        
        function [] = addJob(obj,job)
            obj.jobList{end+1} = job;
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
                
                % setup for output tar via job name
                tmp = strrep(obj.dagline_vars_output_line,'#jobName#',jobName);
                tmp = strrep(tmp,'#N#',num2str(nargs+1)); 
                fprintf(fileID,'%s\n',tmp);
                
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
                
                % setup for untar output in post script
                tmp = strrep(obj.datline_tar_output,'#jobName#',jobName);
                tmp = strrep(tmp,'#outputLocation#',obj.outputLocation);
                fprintf(fileID,'%s\n',tmp);
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
            if ~isempty(obj.prescriptFile)
                scpFileList{6} = [obj.tmpFilesLocation obj.prescriptFile];
            end
        end
        
        function [] = submitDag(obj,varargin)
            if nargin >= 2
                maxidle = varargin{1};
            end
            if nargin >= 3
                maxpost = varargin{2};
            end
            obj.renderDagFile();
            scpList = obj.generate_scpFileList();
            dirCMD_logs_out = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stdout/'''];
            dirCMD_logs_err = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/logs/stderr/'''];
            dirCMD_output = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir -p /home/nate/condorFunctions/#directory#/output/'''];
            [status result] = system(strrep(dirCMD_logs_out,'#directory#',obj.jobFunction));
            [status result] = system(strrep(dirCMD_logs_err,'#directory#',obj.jobFunction));
            [status result] = system(strrep(dirCMD_output,'#directory#',obj.jobFunction));
            dirCMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'mkdir /home/nate/condorFunctions/#directory#/'''];
            [status result] = system(strrep(dirCMD,'#directory#',obj.jobFunction));
            CMD = 'scp -P 50118 #srcfile# nate@128.104.98.118:/home/nate/condorFunctions/#directory#/#desfile#';
            CMD = strrep(CMD,'#directory#',obj.jobFunction);
            for f = 1:numel(scpList)
                [pth nm ext] = fileparts(scpList{f});
                tCMD = strrep(CMD,'#desfile#',[nm ext]);
                tCMD = strrep(tCMD,'#srcfile#',scpList{f});
                [status result] = system(tCMD);
            end

            % submit the job dag
            dagName = obj.generate_dagName();
            CMD = ['ssh -p 50118 nate@128.104.98.118 ''' 'cd /home/nate/condorFunctions/#directory#/; condor_submit_dag -maxidle ' num2str(maxidle) ' -maxpost ' num2str(maxpost) ' ' dagName ''''];
            CMD = strrep(CMD,'#directory#',obj.jobFunction);
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
    end
    
    
    
    
end