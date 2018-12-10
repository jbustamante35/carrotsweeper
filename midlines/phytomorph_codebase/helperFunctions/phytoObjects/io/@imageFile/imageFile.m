classdef imageFile < matlab.mixin.Copyable
    properties
        fileName;
        iticket;
        % make class name 
        %className = 'http://example.org/matlabObjects/imageFile';
    end
    
    properties (Constant)
        % make semantics class name for reference
        %className = [baseDomain.className 'matlabObjects/imageFile'];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = imageFile(varargin)
            obj.fileName = '';
            if nargin == 1
                obj.fileName = varargin{1};
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set file name
        function [] = setFileName(obj,fileName)
            obj.fileName = fileName;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get file name
        function [fileName] = getFileName(obj)
            fileName = obj.fileName;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get file name
        function [fileName] = getFullFileName(obj)
            fileName = getFileName(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read
        function [I] = read(obj)
            I = myReader(obj.fileName);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get file name
        function [nm] = getName(obj)
            [pth,nm,ext] = fileparts(obj.fileName);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % issue ticket for irods
        function [] = issueTicket(obj,uses)
            if nargin == 1
                uses = 1;
            end
            cmd = ['iticket create read ' obj.fileName];
            [o,r] = system(cmd);
            fidx = strfind(r,':');
            obj.iticket = r(fidx+1:end-1);
            cmd = ['iticket mod ' obj.iticket ' uses ' uses];
            [o,r] = system(cmd);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % delete ticket for irods
        function [] = deleteTicket(obj)
            cmd = ['iticket delete ' obj.iticket];
            system(cmd);
        end
        %{
        % removed to redeploy without java
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert to jObject
        function[jObject] = toJobject(obj,store)
            import phytoG.locked.Bpersist.Bfs.abstractions.Buri_file;
            jObject = Buri_file(store);
            %jObject.
        end
        %}
        function [fmoType] = persist(obj,conn)
            %import org.openrdf.model.vocabulary.RDF;
            %import org.openrdf.model.vocabulary.RDFS;
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % create value factory
            vf = conn.getRepository().getValueFactory();
            % make a feature map object
            fmoType = vf.createURI(['irods:/' obj.fileName]);
            % make a class object
            fmoClassName = vf.createURI(obj.className);
            % register the image file
            conn.add(fmoType,RDF.TYPE,fmoClassName,cont);
        end
    end
    
    methods (Static)
        function [] = renderClassType(conn)
            %import org.openrdf.model.vocabulary.RDF;
            %import org.openrdf.model.vocabulary.RDFS;
            %% create value factory
            vf = conn.getRepository().getValueFactory();
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % make grand name space
            nameSpace = 'http://example.org/matlabObjects/';
            % make a feature map object
            fmoClass = vf.createURI(nameSpace,'imageFile');
            conn.add(fmoClass,RDF.TYPE,RDFS.CLASS,cont);
        end
    end
end