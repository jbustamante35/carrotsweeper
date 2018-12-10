classdef iFile < handle & tagable
    
    properties
        remoteFileName;
        localFileName;
    end
    
    properties (Constant)
        localPath = '/mnt/tetra/nate/local_iFile/';
    end
    
    
    methods
        
        function [obj] = iFile(rName)
            [pth,nm,ext] = fileparts(rName);
            obj.localFileName = [iFile.localPath nm ext];
            obj.remoteFileName = rName;
            if obj.existsR()
                obj.loadTags();
            end
        end
        
        function [r] = existsL(obj)
            r = exist(obj.localFileName,'file');
        end
        
        function [r] = existsR(obj)
            CMD = ['ils ' obj.remoteFileName];
            [r,o] = system(CMD);
            if r == 0
                r = true;
            else
                r = false;
            end
        end
        
        function [] = put(obj)
            CMD = ['iput -V ' obj.localFileName ' ' obj.remoteFileName];
            [r,o] = system(CMD);
        end
        
        function [] = get(obj)
            CMD = ['iget -V ' obj.remoteFileName ' ' obj.localFileName];
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