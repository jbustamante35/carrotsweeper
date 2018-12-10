classdef imageStack < myHS_X
    properties        
        iplantUser;
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = imageStack(varargin)
            obj@myHS_X('imageFile');
            if nargin == 1
                if ~isjava(varargin{1})
                    for e = 1:numel(varargin{1})
                        iF = imageFile(varargin{1}{e});
                        putImage(obj,iF);
                    end
                else
                    itr = varargin{1}.iterator();
                    while itr.hasNext()
                        iF = itr.next();
                        iF = char(iF.getFullFileName());
                        putImage(obj,imageFile(iF));
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get image from list
        function [iF] = getImage(obj,n)
            iF = getElement(obj,n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put image into list
        function [] = putImage(obj,iF,n)
            if nargin == 2
                n = numel(obj.S)+1;
            end
            putElement(obj,iF,n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sort
        function [] = sort(obj)
            try
                % loop over all the images
                for e = 1:numel(obj)
                    cur = getElement(obj,e);
                    nm(e) = str2num(cur.getName()); 
                end
                % sort the numbers
                [j sidx] = sort(nm);
                n = obj.copy();                
                for e = 1:numel(obj)
                    cur = getElement(obj,sidx(e));
                    putElement(n,cur,e);
                end
                obj.S = n.S;
            catch ME
                ME;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % size
        function [sz] = size(obj)
            sz = numel(obj.S);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get 
        function [ret] = get(obj,idx)
            ret = obj.S{idx+1};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get
        function [ret] = getLocal(obj,idx)
            ret = get(obj,idx);
            ret = irods2fuse({{ret.getFullFileName()}});            
            ret = imageFile(ret{1}{1});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % xform to iRods syntax
        function [] = xForm_iRods(obj,user,moniker)
            for e = 1:numel(obj.S)
                fn = obj.S{e}.getFullFileName();
                if isFUSE(fn)                     
                    fn = fuse2irods({{fn}},user);
                    fn = fn{1}{1};
                end
                obj.S{e}.setFileName(fn);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % xform to local file
        function [] = xForm_local(obj,path)
            for e = 1:numel(obj)
                [pth nm ext] = fileparts(obj.getImage(e).fileName);
                obj.putImage(imageFile([path nm ext]),e);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % issueTickets over each file in the stack
        function [] = issueTickets(obj,varargin)
            for e = 1:numel(obj)
                tmp = obj.S{e};
                tmp.issueTicket(varargin{:});
                obj.S{e} = tmp;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % persist to irods
        function [] = persist(obj,location)
            %cmd = ['imkdir /iplant/home/' obj.iPlantUser /]
            %save(
        end
        %{
        % removed to deploy on iplant
        function [jObject] = toJobject(obj,store,fileSystem)
            import phytoG.locked.BdataObjects.fileSystem.implementations.imageList;
            import phytoG.locked.Bpersist.Bfs.abstractions.Buri_file;
            jObject = imageList(store);
            jObject.fuse();
            for e = 1:numel(obj)
                tic
                imageFile = Buri_file(store);                
                imageFile.setFileSystem(fileSystem);
                imageFile.setFile(obj.S{e}.fileName);
                %imageFile.persist();
                jObject.insertImage(imageFile);
                toc
            end
        end
        %}
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % import from cell array
    methods (Static)
        function [isS] = imageStackSet(set)
            isS = myHS_X('imageStack');
            for e = 1:numel(set)
                isS{e} = imageStack(set{e});
                isS{e}.sort();
            end
        end
  
    end
end