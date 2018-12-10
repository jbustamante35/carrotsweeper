classdef jobFactory < handle
    properties
        kiosk;
    end
    
    methods
        function [obj] = jobFactory()
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % look for imageList by type
        function [iS] = scanKiosk(obj,type)
            iS = {};
            import phytoG.locked.Bpersist.Bos.implementations.*
            import com.mongodb.*;
            import java.util.Map.*;
            import java.util.HashMap;
            import phytoG.locked.BdataObjects.fileSystem.implementations.imageList;
            % load from object store
            obj.kiosk.accessResource();
            obj.kiosk.setCollection('fileCollection');
            % search map
            qMap = HashMap();
            qMap.put('_Bot','imageList');
            qMap.put('_k_m._pnode._k_m._imageSetType',type);
            cursor = obj.kiosk.search(qMap);
            itr = cursor.iterator();
            % loop over search search results
            c = 1;
            while itr.hasNext()
                imageList = imageList(itr.next(),obj.kiosk);
                if ~logLookUp(obj,type,imageList.getID())
                    iS{c} = imageList;
                    c = c + 1;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % generate jobs
        function [jB] = generateJobs(obj,type)
            iS = scanKiosk(obj,type);
            for e = 1:numel(iS)
                iport = inPort();
                iport.source = obj.kiosk;
                iport.sourceType = 'mongoDB';
                iport.id = char(iS{e}.getID());
                iport.pull();
                oport = outPort();
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        % set the kiosk
        function [] = setKiosk(obj,kiosk) 
            obj.kiosk = kiosk;
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % get the kiosk
        function [k] = getKiosk()
            k = obj.kiosk;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        % clear execute record
        function [] = clearExecuteLog(obj,type)
            executeLog = getLogFromKiosk(obj,type);
            executeLog.clearLog();
            obj.kiosk.accessResource();
            obj.kiosk.setCollection('executeCollection');
            obj.kiosk.put(executeLog);
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % look for execute record via type and id
        % if does not exsist, then create execute record
        function [b] = logLookUp(obj,type,id)
            % import java
            import java.util.HashMap;
            import phytoG.locked.Bcompute.implementations.*;
            % acess the kiosk- mongoDB
            obj.kiosk.accessResource();
            obj.kiosk.setCollection('executeCollection');
            % setup search and perform
            qMap = HashMap();
            qMap.put('_Bot','Bexecute_log');
            qMap.put('_k_m._pnode._k_m._methodName',type);
            qMap.put('_k_m._dnode._k_l',id);
            cursor = obj.kiosk.search(qMap);
            itr = cursor.iterator();
            if ~(itr.hasNext())
                b = 0;
            else
                b = 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % look for execute record
        function [b] = logExist(obj,type)
            % import java
            import java.util.HashMap;
            import phytoG.locked.Bcompute.implementations.*;
            % acess the kiosk- mongoDB
            obj.kiosk.accessResource();
            obj.kiosk.setCollection('executeCollection');
            % setup search and perform
            qMap = HashMap();
            qMap.put('_Bot','Bexecute_log');
            qMap.put('_k_m._pnode._k_m._methodName',type);
            cursor = obj.kiosk.search(qMap);
            itr = cursor.iterator();
            b = itr.hasNext();
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % delete execute log
        function [] = deleteExecuteLog(obj,type)
            if logExist(obj,type)
                ret = getLogFromKiosk(obj,type);
                obj.kiosk.delete(ret.getID());
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % get execution record
        function [ret] = getLogFromKiosk(obj,type)
             % import java
            import java.util.HashMap;
            import phytoG.locked.Bcompute.implementations.*;
            % acess the kiosk- mongoDB
            obj.kiosk.accessResource();
            obj.kiosk.setCollection('executeCollection');
            % setup search and perform
            qMap = HashMap();
            qMap.put('_Bot','Bexecute_log');
            qMap.put('_k_m._pnode._k_m._methodName',type);
            cursor = obj.kiosk.search(qMap);
            itr = cursor.iterator();
            if (itr.hasNext())
                ret = itr.next();
                ret = Bexecute_log(ret,obj.kiosk);
            else
                ret = createExecuteLog(obj,type);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % post execute log
        function [] = postLog(obj,exeLog)
            obj.kiosk.put(exeLog);
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        % create execute record
        function [eLog] = createExecuteLog(obj,type)
            import phytoG.locked.Bcompute.implementations.*;
            b = logExist(obj,type);
            if ~b
                eLog = Bexecute_log();
                eLog.setMethodName(type);
                obj.kiosk.put(eLog);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%{

%}