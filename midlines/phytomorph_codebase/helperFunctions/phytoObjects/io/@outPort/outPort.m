classdef outPort < matlab.mixin.Copyable
    properties
        % authentication
        userName;
        % for file system io
        writeBasePath;
        uniqueKey;
        specialCharater = '--';
        nDepth = 5;
    end
    methods
        
        function [obj] = outPort()
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the write path base
        function [] = setWriteBase(obj,writebase)
            obj.writeBasePath = writebase;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the write path base
        function [writebase] = getWriteBase(obj)
             writebase = obj.writeBasePath;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set user for outPath
        function [] = setUser(obj,user)
            obj.userName = user;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate unique key
        function [uK] = genUniqueKey(obj,filebase)
            obj.uniqueKey = genUQkey(filebase,obj.specialCharater,obj.nDepth);
            obj.writeBasePath = [obj.writeBasePath obj.uniqueKey filesep];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push
        function [fileList] = push(obj,fileName,data)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % augment the path with the name ext
            [spath,nm,ext] = fileparts(fileName);
            if ~isempty(spath)
                kPath = [obj.writeBasePath spath filesep];
            else
                kPath = [obj.writeBasePath];
            end
            
            
            obj.generatePath(kPath)
            kpPath = [kPath nm ext];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if isFUSE
            if isFUSE(obj.writeBasePath)
                % transform
                fileList = fuse2irods({{kpPath}},obj.userName);
                kpPath = fileList{1}{1};
                OLDfile = kpPath;
                [spath,nm,ext] = fileparts(kpPath);
                kpPath = ['/tmp' filesep nm ext];
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % assign for writer
            out.fileName = kpPath;
            out.d = data;
            % call to write data
            fileList = myWriter(kpPath,out);
            fileList = fileList{1};
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if isFUSE
            if isFUSE(obj.writeBasePath)                
                xfer_put(kpPath,OLDfile);
                delete(kpPath);
                fileList = OLDfile;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create base path
        function [] = generatePath(obj,genPath)
            if isFUSE(genPath)
                % transform
                fileList = fuse2irods({{genPath}},obj.userName);
                genPath = fileList{1}{1};
            end
            if isIRODS(genPath)
                cmd = ['imkdir -p "' genPath '"'];
                [r,o] = system(cmd);
            else
                mkdir(genPath);
            end
        end
        
    end
    
end
