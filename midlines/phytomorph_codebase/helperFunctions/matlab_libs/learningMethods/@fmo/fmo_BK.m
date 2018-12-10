classdef fmo < handle
    % feature map object: fmo
    % a matrix form of the feature information
    % vector if raw image
    % list of row vectors else
    properties
        % the point detector object : featureMapBank+featureFunctionBank
        ptde = [];
        % the feature map data for this object
        featureMap = [];
        % the index map to the domain
        indexMap = [];
        % the size of the feature map
        rawSize = [];
        % the size of the patches
        featurePatchSize = [];
        % the human readable moniker for the feature map
        moniker = [];
        % the unique key for the featuremap
        key = [];
        % the type key for the feature map
        typeKey = '';
        % head node name : the unique file name for the head node
        headNodeName = '';
        
        % make semantics class name for reference
        className = 'http://example.org/matlabObjects/featureMapObject';
    end
    
    methods
        function [obj] = fmo(varargin)
            if nargin >= 1
                obj.ptde = varargin{1};
            end
            if nargin >= 2
                obj.moniker = varargin{2};
            end
            % if this then register the type : needed
            if nargin >= 3
                %old random key
                %obj.typeKey = ['t_' num2str(string2hash([datestr(clock) '-' num2str(randi(10000,1))]))];
                % new non-random key
                obj.typeKey = ['t_' num2str(string2hash(obj.moniker))];
                obj.ptde.registerFeatureMapType(obj);
            end
            % delay the generation of the unique key until persistance time
            % or until deposit time
            % obj.key = obj.generateObjectKey();
            
            
            % make grand name space
            %nameSpace = 'http://example.org/productionSystem/';
            % make a feature map object
            %obj.className = vf.createURI(nameSpace,'featureMapObject');
        end
        
        function [] = setHeadNodeName(obj,headNodeName)
            obj.headNodeName = headNodeName;
        end
        
        function [headNodeName] = getHeadNodeName(obj)
            if ~isempty(obj.headNodeName)
                headNodeName = obj.headNodeName;
            else
                sPath = obj.ptde.oPath;
                load([sPath obj.key '.mat'],'headNodeName');
            end
        end
        
        function [nKey] = generateObjectKey(obj)
            % random+clock key
            toHash = [datestr(clock) '-' num2str(randi(10000,1))];
            nKey = ['o_' num2str(string2hash(toHash))];
            fprintf(['Hashed:' toHash '-->' nKey '\n'])
        end
        
        function [] = setptde(obj,ptde)
            obj.ptde = ptde;
        end
        
        function [] = setMoniker(obj,moniker)
            obj.moniker = moniker;
        end
        
        function [] = setData(obj,data)
            obj.featureMap = data;
        end
        
        function [featureMap] = getData(obj,index)
            tm = clock;
            fprintf(['start loading call for key @' obj.key '\n']);
            sPath = obj.ptde.oPath;
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
            fprintf(['end loading call for key @' obj.key ':' num2str(etime(clock,tm)) '\n']);
        end
        
        function [SZ] = getFeatureMapSize(obj)
             sPath = obj.ptde.oPath;
             datFile = [sPath obj.key '.h5'];
             SZ = h5readatt(datFile,'/data','size');
        end
        
        function [] = setIndexMap(obj,index)
            obj.indexMap = index;
        end
        
        function [indexMap] = getIndexMap(obj)
            sPath = obj.ptde.oPath;
            load([sPath obj.key '.mat'],'indexMap');
        end
        
        function [] = setRawSize(obj,sz)
            obj.rawSize = sz;
        end
        
        function [rawSize] = getRawSize(obj)
            sPath = obj.ptde.oPath;
            load([sPath obj.key '.mat'],'rawSize');
        end
        
        function [] = setPatchSize(obj,psz)
            obj.featurePatchSize = psz;
        end
        
        function [featurePatchSize] = getPatchSize(obj)
            sPath = obj.ptde.oPath;
            load([sPath obj.key '.mat'],'featurePatchSize');
        end
        
        function [] = setTypeKey(obj,typeKey)
            obj.typeKey = typeKey;
        end
        
        function [typeKey] = getTypeKey(obj)
            if isempty(obj.typeKey)
                sPath = obj.ptde.oPath;
                load([sPath obj.key '.mat'],'typeKey');
            else
                typeKey = obj.typeKey;
            end
        end
        
        function [] = setKey(obj)
            obj.key = ['hnn_' obj.getHeadNodeName '_typeKey_' obj.getTypeKey()];
        end
        
        function [] = persist(obj)
            obj.setKey();
            % get the save path for the object from the point detector
            % object
            sPath = obj.ptde.oPath;
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
            
            
            rawSize = obj.rawSize;
            obj.rawSize = [];
            
            
            featurePatchSize = obj.featurePatchSize;
            obj.featurePatchSize = [];
            
            
            moniker = obj.moniker;
            obj.moniker = '';
            
            headNodeName = obj.headNodeName;
            obj.headNodeName = '';
            
            typeKey = obj.typeKey;
            obj.typeKey = '';
            
            
            % get the key for saving
            key = obj.key;
            
            fprintf(['start persisting@' moniker '@' obj.key '\n']);tic
            matFile = [sPath obj.key '.mat'];
            datFile = [sPath obj.key '.h5'];
            save(matFile,'indexMap','rawSize','featurePatchSize','moniker','key','typeKey','headNodeName');
            h5create(datFile,'/data',size(featureMap));
            h5write(datFile,'/data',featureMap);
            h5writeatt(datFile,'/data','size',size(featureMap));
            fprintf(['end persisting@' moniker '@' obj.key ':' num2str(toc) '\n']);
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
            % make grand name space
            psNameSpace = 'http://example.org/productionSystem/';
            % make a feature map object
            fmoClass = vf.createURI(psNameSpace,'featureMapObject');
            conn.add(fmoClass,RDF.TYPE,RDFS.CLASS,cont);
        end
    end
end