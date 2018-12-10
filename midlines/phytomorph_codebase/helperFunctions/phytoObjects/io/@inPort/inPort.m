classdef inPort < matlab.mixin.Copyable
    properties
        loginName;
        fileList;
        % for pull
        source;
        sourceType;
    end
    methods
        
        function [obj] = inPort()
            
        end
        
        function [] = issueTickets(obj,varargin)
            obj.fileList.issueTickets(varargin{:});
        end
        
        function [] = setLogin(obj,loginName)
        
        end
        
        function [] = getLogin(obj,loginName)
        
        end
        
        function [] = localizeFileList(obj,path)
            if nargin == 1
                path = './';
            end
            obj.fileList.xForm_local(path);
        end
        
        function [] = pull(obj)
            switch obj.sourceType
                % obj.source = path to scan
                case 'structuredPathScan'
                    ret = structuredPathScan(obj.source);
                    if numel(ret) ~= 0
                        obj.fileList = ret{1};
                    else
                        obj.fileList = ret;
                    end
                    
                % obj.source is list of files
                case 'fileList'
                    obj.fileList = imageStack(obj.source);
                % obj.source is mongoDB connection
                % id fileList to load
                %{
                case 'mongoDB'
                    import java.util.HashMap;
                    import phytoG.locked.BdataObjects.fileSystem.implementations.imageList;
                    obj.source.db.accessResource();
                    obj.source.db.setCollection('fileCollection');
                    qMap = HashMap();
                    qMap.put('_id',obj.source.id);
                    cursor = obj.source.search(qMap);
                    itr = cursor.iterator();
                    iL = itr.next();
                    iL = imageList(iL,obj.source);
                    iL = iL.getFileList();
                    obj.fileList = imageStack(iL);
                %}
            end
        end
        
    end
end