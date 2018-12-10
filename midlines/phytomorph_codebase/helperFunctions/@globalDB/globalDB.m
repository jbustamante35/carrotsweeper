classdef globalDB < handle & tagable
    
    properties
        id;
        isMEM;
        dbLocalFile;
        dbRemoteFile;
    end
    
    
    properties (Constant)
        localPath = '/mnt/scratch4/localDB/';
        remotePath = '/iplant/home/nmiller/publicData/globalDB/';
    end
    
    
    
    methods
        
        
        function [obj] = globalDB(varargin)
            % if no file name is given-then default to in memory database
            if nargin == 0
                obj.dbFile = ':memory:';
                obj.isMEM = true;
            else
                dbFile = varargin{1};
                tmpLocal = globalDB.localPath;
                tmpRemote = globalDB.remotePath;
                tmpLoad = 'Q';
                if nargin >= 2
                    if ~isempty(varargin{2})
                        tmpLocal = varargin{2};
                    else
                        tmpLocal = globalDB.localPath;
                    end
                end
                if nargin >= 3
                    if ~isempty(varargin{3})
                        tmpRemote = varargin{3};
                    else
                        tmpRemote = globalDB.remotePath;
                    end
                end
                if nargin >= 4
                    toLoad = varargin{4};
                end
                % is MEM = false
                obj.isMEM = false;
                obj.dbLocalFile = [tmpLocal dbFile];
                obj.dbRemoteFile = [tmpRemote dbFile];
                obj.dbFile = obj.dbLocalFile;
                % 
                if obj.existsR()
                    switch toLoad
                        case 'Q'
                            
                        case 'N'
                            
                        case 'Y'
                    end
                    fprintf(['dB exists@remote.\n']);
                end
            end
            % opeN
            mksqlite('open',obj.dbFile);
            % set blob type
            mksqlite('typedBLOBs', 1 );
            % wrapping
            mksqlite('param_wrapping', 1 );
        end
        
        
        
        function [] = publish(obj)
            if ~strcmp(obj.dbFile,':memory:')
                cmd = ['iput -V ' obj.dbFile ' ' obj.dbRemoteFile];
            end
        end
        
        
        function [r] = existsR(obj)
            CMD = ['ils ' obj.dbRemoteFile];
            [r,o] = system(CMD);
            if r == 0
                r = true;
            else
                r = false;
            end
        end
        
        
        function [] = iput(obj)
            mksqlite('close');
            CMD = ['iput -V ' obj.dbLocalFile ' ' obj.dbRemoteFile];
            [r,o] = system(CMD);
            mksqlite('open',obj.dbLocalFile);
        end
        
        function [] = iget(obj)
            CMD = ['iget -V ' obj.dbRemoteFile ' ' obj.dbLocalFile];
            [r,o] = system(CMD);
        end
        
        function [] = loadTags(obj)
            CMD = ['imeta ls -d ' obj.remoteFileName];
            [r,o] = system(CMD);
            tidx = strfind(o,char(10));
            o(1:(tidx(1)-1)) = [];
            o = ['----' o(1:end-1) o(end) '----' o(end)];
            tidx = strfind(o,['----' char(10)]);
            for e = 1:(numel(tidx)-1)
                chop = o((tidx(e)+5):tidx(e+1)-2);
                k = strfind(chop,'attribute:');
                v = strfind(chop,'value:');
                u = strfind(chop,'units:');
                k = chop(12:(v(1)-2));
                v = chop((v(1)+7):(u(1)-2));
                obj.setTag(k,v,false);
            end
        end
        
        
        
        
        
    end
    
    
    
end