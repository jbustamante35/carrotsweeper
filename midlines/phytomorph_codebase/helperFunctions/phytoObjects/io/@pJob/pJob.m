classdef pJob < matlab.mixin.Copyable
    properties
        % job type
        methodType;
        
        % inPorts for read and write 
        inPort;
        outPort;
        
        % functions which need to be run
        Sfunc;
        func;
        para;
        
        % display during job
        disp;
        
        % ticket
        iticket;
        
        % results container
        results;
    end
    
    methods
        function [obj] = pJob()
            obj.disp = 0;
        end
        
        function [] = setIport(obj,inPort)
            obj.inPort = inPort;
        end
        
        function [] = setOport(obj,outPort)
            obj.outPort = outPort;
        end
        
        function [] = setResultsContainer(obj,rc)
            obj.results = rc;
        end
            
        function [] = setDisp(obj,disp)
            obj.disp = disp;
        end
      
        function [] = attachParaField(obj,key,value)
            obj.para = setfield(obj.para,key,value);
        end
        
        function [] = attachFunction(obj,func)
            obj.func(end+1).f = func;
        end
        
        function [] = attachStagingFunction(obj,sFunc)
            obj.Sfunc(end+1).f = sFunc;
        end
        
        function [] = executeStagingFunction(obj)
            % call each function in the job package
            for fn = 1:numel(obj.Sfunc)
                out(fn).val = obj.Sfunc(fn).f(obj);
            end
        end
        
        function [] = issueTickets(obj,varargin)
            obj.inPort.issueTickets(varargin{:});
        end
        
        function [out] = run(obj)
            % call each function in the job package
            for fn = 1:numel(obj.func)
                try
                    out(fn).val = obj.func(fn).f(obj);
                catch ME
                    getReport(ME)
                end
            end
        end
        
        function [] = postJob(obj,varargin)
            % create local file name
            fn = [tempname '.mat'];
            % save the file locally
            save(fn,'obj');
            % create the remote path
            postfn = ['/iplant/home/' varargin{1} '/phytoMorph/jobs/'];
            % make directory for remote path
            cmd = ['imkdir ' postfn];
            system(cmd);
            [pth nm ext] = fileparts(fn);
            postfn = [postfn nm ext];
            % transfer the job file to the remote path
            cmd = ['iput ' fn ' ' postfn];
            system(cmd);
            
            % issue once only read ticket for the job
            cmd = ['iticket create read ' postfn];
            [o,r] = system(cmd);
            fidx = strfind(r,':');
            obj.iticket = r(fidx+1:end-1);
            cmd = ['iticket mod ' obj.iticket ' uses 1'];
            [o,r] = system(cmd);
            
            cmd = ['imeta add -d ' postfn ' program atmoWorker'];
            [o,r] = system(cmd);
            cmd = ['imeta add -d ' postfn ' status pending'];
            [o,r] = system(cmd);
            
        end
          
        function [] = generateCondorSubmit(obj,packFile)
            fileID = fopen('/mnt/scratch1/junk.submit','w');
            l{1} = 'universe = vanilla \n';
            l{2} = 'should_transfer_files = YES \n';
            l{3} = 'transfer_input_files = ';
            l{4} = 'Executable = condor_worker\n';
            l{5} = 'requirements = (Arch == "x86_64") && (OpSys == "LINUX") && (Disk > 500000) && (TotalMemory >= (3*1024))\n';
            l{6} = '+AccountingGroup = "spalding\n';
            l{7} = '+Group = "spalding"\n';           
            l{8} = 'priority = 13\n';
            l{9} = 'output = ';
            l{10}= 'error = ./algo.error\n';
            l{11}= 'Queue\n';
            
            INLINE = 3;
            for e = 1:numel(obj.inPort.fileList)
            
            end
            
            for e = 1:numel(l)
                fprintf(fileID,l{i},'char');
            end
            
            fclose(fileID);
            obj.inPort.localizeFileList();
        end
         
    end
end