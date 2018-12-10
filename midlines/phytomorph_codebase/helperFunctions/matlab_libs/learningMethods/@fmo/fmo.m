classdef fmo < handle
    % feature map object: fmo
    % a matrix form of the feature information
    % vector if raw image
    % list of row vectors else
    properties
        % the feature map data for this object
        featureMap = [];
        % the index map to the domain
        indexMap = [];
        % the size of the feature map
        rawSize = [];
        % the size of the patches
        featurePatchSize = [];
        % the human readable moniker for the feature map
        moniker = '';
        % the unique key for the featuremap
        key = [];
        
        % data location
        oPath = '';
    end
    
    properties (Constant)
        % make semantics class name for reference
        className = [baseDomain.className 'matlabObjects/featureMapObject'];
        subDomain = 'matlabObjects/';
        % make representable tag
        representable = [baseDomain.className 'matlabObjects/representable'];
    end
    
    methods
        function [obj] = fmo(varargin)
            % delay the generation of the unique key until persistance time
            % or until deposit time
            obj.key = obj.generateObjectKey();
        end
        
        function [nKey] = generateObjectKey(obj)
            % random+clock key
            toHash = [datestr(clock) '-' num2str(randi(10000,1))];
            nKey = ['o_' num2str(string2hash(toHash))];
            fprintf(['Hashed for feature map object:' toHash '-->' nKey '\n'])
        end
        
        function [] = setMoniker(obj,moniker)
            obj.moniker = moniker;
        end
        
        function [] = setData(obj,data)
            obj.featureMap = data;
        end
        
        function [featureMap] = getData(obj,index)
            tm = clock;
            fprintf(['fmo: start loading call for key @' obj.key '\n']);
            sPath = obj.getoPath();
            datFile = [sPath obj.key '.h5'];
            if nargin > 1
                if ~iscell(index)
                    im = obj.getIndexMap();
                    [~,sidx,~] = intersect(im,index);
                    featureMap = h5read(datFile,'/data');
                    featureMap = featureMap(:,sidx);
                else
                    SZ = h5readatt(datFile,'/data','size');
                    STR = [(SZ(1)-1)/2+1 1];
                    OFF = [1 SZ(2)];
                    featureMap = h5read(datFile,'/data/',STR,OFF);
                end
            else
                featureMap = h5read(datFile,'/data');
            end
            fprintf(['fmo: end loading call for key @' obj.key ':' num2str(etime(clock,tm)) '\n']);
        end
        
        function [SZ] = getFeatureMapSize(obj)
             sPath = obj.getoPath();
             datFile = [sPath obj.key '.h5'];
             SZ = h5readatt(datFile,'/data','size');
        end
        
        function [] = setIndexMap(obj,index)
            obj.indexMap = index;
        end
        
        function [indexMap] = getIndexMap(obj)
            sPath = obj.getoPath();
            load([sPath obj.key '.mat'],'indexMap');
        end
        
        function [] = setRawSize(obj,sz)
            obj.rawSize = sz;
        end
        
        function [rawSize] = getRawSize(obj)
            sPath = obj.getoPath();
            load([sPath obj.key '.mat'],'rawSize');
        end
        
        function [] = setPatchSize(obj,psz)
            obj.featurePatchSize = psz;
        end
        
        function [featurePatchSize] = getPatchSize(obj)
            sPath = obj.getoPath();
            load([sPath obj.key '.mat'],'featurePatchSize');
        end
        
        function [] = setKey(obj)
            obj.key = ['hnn_' obj.getHeadNodeName '_typeKey_' obj.getTypeKey()];
        end
        
        function [oPath] = getoPath(obj)
            oPath = obj.oPath;
        end
        
        function [fileObject] = persist(obj,pInfo)
            obj.oPath = pInfo.oPath;
            fileObject = obj.registerData(pInfo);
            obj.persistData(pInfo);
        end
        
        function [] = persistData(obj,pInfo)
            
            % get the save path for the object data
            sPath = pInfo.oPath;
            %%%%%%%%%%%
            %feature map
            %%%%%%%%%%%
            % get the feature map for the object
            featureMap = obj.featureMap;
            % clear the feature map
            obj.featureMap = [];
            
            % get the index map in savable form
            indexMap = obj.indexMap;
            obj.indexMap = [];
            
            % get the raw size
            rawSize = obj.rawSize;
            obj.rawSize = [];
            
            % get the feature map size
            featurePatchSize = obj.featurePatchSize;
            obj.featurePatchSize = [];
            
            % get the moniker
            moniker = obj.moniker;
            obj.moniker = '';
            
            % get the key for saving
            key = obj.key;
            uniqueKey = num2str(string2hash(obj.key));
            
            
            fprintf(['start persisting@' moniker '@' obj.key '\n']);tic
            matFile = [sPath obj.key '.mat'];
            datFile = [sPath obj.key '.h5'];
            save(matFile,'indexMap','rawSize','featurePatchSize','moniker','key','obj');
            %save(matFile,'obj');
            h5create(datFile,'/data',size(featureMap));
            h5write(datFile,'/data',featureMap);
            h5writeatt(datFile,'/data','size',size(featureMap));
            fprintf(['end persisting@' moniker '@' obj.key ':' num2str(toc) '\n']);
        end
        
        function [fileObject] = registerData(obj,pInfo)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            oPath = pInfo.oPath;
            conn = pInfo.conn;
            vf = conn.getRepository().getValueFactory();
            
            % spool the object
            %uniqueKey = num2str(string2hash(obj.key));
            uniqueKey = obj.key;
            
            % register the object
            fileNameSpace = ['file:/' oPath];
            % make a feature map object
            fileObject = vf.createURI(fileNameSpace,[uniqueKey '.mat']);
        end
        
        function [] = clear(obj)
          
        end
        
        function [] = view(obj)
            tmpD = obj.getData({});
            viewIndex = (size(tmpD,1)-1)/2 + 1;
            tmpD = tmpD(viewIndex,:);
            tmpD = reshape(tmpD,obj.getRawSize());
            imshow(tmpD,[]);
        end
    end
    
    methods (Static)
        function [] = renderClassType(conn)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            %% create value factory
            vf = conn.getRepository().getValueFactory();
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % make a feature map object
            fmoClass = vf.createURI(fmo.className);
            conn.add(fmoClass,RDF.TYPE,RDFS.CLASS,cont);
        end
        
        function [fmoSubClass] = renderSubClassType(conn,subType,toRepresent)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            % create value factory
            vf = conn.getRepository().getValueFactory();
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % make a feature map object classname
            fmoClass = vf.createURI(fmo.className);
            % make a feature map object
            fmoSubClass = vf.createURI([baseDomain.className fmo.subDomain],subType);
            % make sub class
            conn.add(fmoSubClass,RDFS.SUBCLASSOF,fmoClass,cont);
            if toRepresent
                % make feature map object representable
                representable = vf.createURI(fmo.representable);
                conn.add(fmoSubClass,RDF.PROPERTY,representable,cont);
            end
        end
        
        function [] = renderSubClassObject(conn,subType,object)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
             % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % create value factory
            vf = conn.getRepository().getValueFactory();
            % make a feature map object
            fmoSubClass = vf.createURI([baseDomain.className fmo.subDomain],subType);
            % render the subClass type PRE RENDER THE TYPE ON PERSISTANCE
            %fmoSubClass = fmo.renderSubClassType(conn,subType);
            % render object as subTypeClass type
            conn.add(object,RDF.TYPE,fmoSubClass,cont);
          
        end
        
        function [type] = getSubClassObjectByType(conn,type)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            vf = conn.getRepository().getValueFactory();
            % make a feature map object
            type = vf.createURI([baseDomain.className fmo.subDomain],type);
        end
        
        function [obj] = loadFromURI(uriList)
            obj = {};
            if ~iscell(uriList)
                tmp = uriList;
                uriList = {};
                uriList{1} = tmp;
            end
            for e = 1:numel(uriList)
                fileName = char(uriList{e}.toString());
                fileName = fileName(7:end);
                tmp = load(fileName,'obj');
                obj{end+1} = tmp.obj;
            end
        end
    end
end