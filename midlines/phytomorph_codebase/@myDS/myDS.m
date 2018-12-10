classdef myDS < handle & matlab.io.Datastore & matlab.io.datastore.Partitionable & globalDB
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % datastore: holds nodes and edges
    % this is starting to rebuild the datastore that I need
    % the basis particle is a node/edge object - enode or nedge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % each basic particle (en:edge-node) has potential for
    %       1: source, target, and data.
    %           --also will store date and key for each en
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % there are four fundamental particles
    % when these particles are all attached together they create a macroparticle complex
    % which will be called tau.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % nu: this particle holds the lambda, pi and delta particles together
    % as such there are many ways to attach the complex - here as of June, 28 2018
    % i will have each particle in the complex point to the nu-center
    %%%%%%%
    % lambda: this particle is responsible for linking the complex to other macro-particle
    % complex(s).  the lambda is attached to a eLambda and iLambda complex
    % the eLambda handles the efflux attachments to other nu particles
    % the iLambda handles the infllux attachments to other nu particles
    %%%%%%%
    % delta: this particle will house the data. as simple as that. for now
    %%%%%%%
    % pi: this particle handles the properties for the data/associated node
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % link types are direct and indirect - indirect are matter and
    % non-matter - non-matter is 
    
    
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % main data store objects
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dbFile
        dbFile;
        % data store
        ds;
        % current read location
        readPtr;
        % current insert location
        emptyPtr;
        
        
         
        adjTables = {};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sequence path through the tree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sequence
        sequence;
        % sequence pointer
        seqPtr = 1;
        % readThrough flag
        toReadThrough = true;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read/load file/image data
        loadFlag;
        base_outPath = '';
        
        
        filesPerPartitions = 10;
        
        effluxRender = true;
        influxRender = true;
        
    end
    
    properties (Constant)
        IMAGE_EXT_LIST = {'tiff','jpg','jpeg','png','bmp','nd2','czi','nms','nef'};
        KEYWORD = {'corn','maize','ear','kernel','cob','carrot','stomata','arabidopsis','seed',...
                    'gravitropism','root','seedling','hypocotyl','potato','mutatnt','cvi'};
        META_EXT_LIST = {'json'};
        BASE = 2.^(0:7);
        
        nonDuplicate = true;
        TYPE = '_type';
        
        NODE = '_node';
        LINK = '_link';
        
        % NOT NEEDED?
        TAU = '_ta';
        % center fu particle
        NU = '_nu';
        % lambda (link particle)
        LAMBDA = '_la';
        % pi (properties particle) 
        PI = '_pi';
        % delta (data particle)
        DELTA = '_de';
        
        
        % matter non-direct link
        NDML = '_ndml';
        
        % efflux and influx node types
        inFluxNode = '_if';
        efFluxNode = '_ef';
        
        % for setting properties
        FLAGKEY = '_fk';
        FLAGVALUE = '_vk';
        
        
        
        
        
        
        
        QUEUE = '_queue';
        QUEUEELE = '_element';
        CONTEXT = '_context';
        
       
        
    end

    methods 

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % basic: init,put(s),get(s),modify(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = myDS(varargin)
            obj@globalDB(varargin{:});
            %obj.adjTables{1} = sparse([],[],[],2^32,2^32,1000);
            %obj.adjTables{2} = sparse([],[],[],2^32,2^32,1000);
            %obj.adjTables{3} = sparse([],[],[],2^32,2^32,1000);
            obj.adjTables{1} = sparse(1000,1000);
            obj.adjTables{2} = sparse(1000,1000);
            obj.adjTables{3} = sparse(1000,1000);
            
            %{
            % if no file name is given-then default to in memory database
            if nargin == 0
                obj.dbFile = ':memory:';
                % open database
                mksqlite('open',obj.dbFile);
            end
            % 
            if ~strcmp(varargin{1},':memory:')
                obj.dbFile = [varargin{1}];
                mksqlite('open',obj.dbFile);
            end
           
            % set blob type
            mksqlite( 'typedBLOBs', 1);
            %}
            
            % create flex-links table
            mksqlite('CREATE TABLE flinks (id INTEGER PRIMARY KEY AUTOINCREMENT,key,time,data,type,target_key,source_key)');
           
            mksqlite('CREATE TABLE tlinks (id INTEGER PRIMARY KEY AUTOINCREMENT,key,time,linkType,parent)');
            
            
            
            
            here = 1;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % general search for rows - return key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [q] = search(obj,fld,con)
            sql = ['SELECT key FROM flinks WHERE '];
            for e = 1:numel(fld)
                if ischar(con{e})
                    sql = [sql fld{e} '=''' con{e} ''' '];
                else
                    sql = [sql fld{e} '=' con{e} ' '];
                end
                sql = [sql 'AND '];
            end
            sql(end-4:end) = [];
            fprintf(['issuing SQL:' sql '\n']);
            q = mksqlite(sql);
        end
        
        function [N] = numberEntries(obj)
            q = mksqlite('SELECT COUNT(*) from flinks');
            N = q.COUNT___;
        end
        
        function [keyPath] = searchTypeNetwork(obj,queryType)
            queryType = lower(queryType);
            if ~strcmp(queryType(1),'_')
                queryType = ['_' queryType];
            end
            q = obj.search({'type'},{[queryType myDS.NU myDS.NODE]});
            sourceKey = q.key;
            N = obj.numberEntries();
            targetVector = 1:N;
            keyPath = obj.dijkstra(sourceKey,targetVector,2);
        end
        
        function [keyPath] = dijkstra(obj,sK,tK,N)
            keyPath = {};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~iscell(sK) && ~isnumeric(sK)
                sK = {sK};
            end
            if ~isnumeric(sK)
                for e = 1:numel(sK)
                    sR(e) = obj.getRowid(sK{e});
                end
            else
                sR = sK; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~iscell(tK) && ~isnumeric(tK)
                tK = {tK};
            end
            if ~isnumeric(tK)
                for e = 1:numel(tK)
                    tR(e) = obj.getRowid(tK{e});
                end
            else
                tR = tK; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            cnt = 1;
            for e1 = 1:numel(sR)
                for e2 = 1:numel(tR)

                    path = dijkstra(obj.adjTables{N},tR(e2),sR(e1));

                    if ~isempty(path)
                        for p = 1:numel(path)
                            keyPath{cnt}{p} = obj.getKey(path(p));
                        end
                        cnt = cnt + 1;
                    end
                end
            end
            
           
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % parameterized get sub particle for node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [sub] = getXNode(obj,key,X)
            switch X
                case 'ef'
                    sub = obj.getEffluxNode(key);
                case 'in'
                    sub = obj.getInfluxNode(key);
                case 'pi'
                    sub = obj.getPiNode(key);
                case 'delta'
                    sub = obj.getDeltaNode(key);
                case 'lambda'
                    
                case 'nu'
                    sub = key;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get efFlux object for node keyList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [efFluxKey] = getEffluxNode(obj,key)
            lam_key = obj.search({'target_key','type'},{key,[myDS.LAMBDA myDS.NODE]});
            efFluxKey = obj.search({'target_key','type'},{lam_key.key,[myDS.efFluxNode myDS.NODE]});
            efFluxKey = efFluxKey.key;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get inFlux object for node keyList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [inFluxKey] = getInfluxNode(obj,key)
           lam_key = obj.search({'target_key','type'},{key,[myDS.LAMBDA myDS.NODE]});
           inFluxKey = obj.search({'target_key','type'},{lam_key.key,[myDS.inFluxNode myDS.NODE]});
           inFluxKey = inFluxKey.key;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get delta node for key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [deltaKey] = getDeltaNode(obj,key)
            deltaKey = obj.search({'target_key','type'},{key,[myDS.DELTA myDS.NODE]});
            deltaKey = deltaKey.key;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get pi node for key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [piKey] = getPiNode(obj,key)
            piKey = obj.search({'target_key','type'},{key,[myDS.PI myDS.NODE]});
            piKey = piKey.key;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get lambda node for key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [lambdaKey] = getLambdaNode(obj,key)
            lambdaKey = obj.search({'target_key','type'},{key,[myDS.LAMBDA myDS.NODE]});
            lambdaKey = lambdaKey.key;
        end
        
        
        function [pointerStruct] = getInfluxLinkList(obj,nuKey)
            inN = obj.getInfluxNode(nuKey);
            pointerStruct = obj.search({'target_key','type'},{inN,myDS.LINK});
            for p = 1:numel(pointerStruct)
                % crawl the link to the object pointed at
                tmpKey = obj.linkCrawlReverse(pointerStruct(p).key,1);
                % get the type of the source object
                tmpType = obj.getType(tmpKey{end});
                % source key
                tmpSourceKey = tmpKey{end};
                % if object entering influx particle is pointing from efflux node
                if strcmp(tmpType,[myDS.efFluxNode myDS.NODE])
                    % forward crawl to center
                    tmpKey = obj.linkCrawl(tmpKey{end},2);
                    tmpKey = tmpKey{end};
                    % get center type
                    tmpType = obj.getType(tmpKey);
                    % default to non-matter - assume link crawl is done
                    tmpLinkType = 'nonmatter';
                    % target key
                    tmpSourceKey = tmpKey;
                    % if type is NDML
                    if strcmp(tmpType,[myDS.NDML myDS.NU myDS.NODE])
                        % get the efflux node of NDML
                        tmpKey = obj.getInfluxNode(tmpKey);
                        % get the pointer entering NDML
                        tmpKey = obj.search({'target_key','type'},{tmpKey,myDS.LINK});
                        tmpKey = tmpKey.key;
                        % crawl to its target - assume isBound
                        tmpKey = obj.linkCrawlReverse(tmpKey,1);
                        tmpKey = obj.linkCrawl(tmpKey{end},2);
                        tmpKey = tmpKey{end};
                        % 
                        tmpType = obj.getType(tmpKey);
                        tmpLinkType = 'matter';
                        % target key
                        tmpSourceKey = tmpKey;
                    end
                end
                pointerStruct(p).targetKey = nuKey;
                pointerStruct(p).type = tmpType;
                pointerStruct(p).linkType = tmpLinkType;
                pointerStruct(p).sourceKey = tmpSourceKey;
            end
            pointerStruct = orderfields(pointerStruct,[5 1 2 4 3]);
        end
        
        
        function [pointerStruct] = getEffluxLinkList(obj,nuKey)
            efN = obj.getEffluxNode(nuKey);
            pointerStruct = obj.search({'source_key','type'},{efN,myDS.LINK});
            for p = 1:numel(pointerStruct)
                % crawl the link to the object pointed at
                tmpKey = obj.linkCrawl(pointerStruct(p).key,1);
                % get the type of the target object
                tmpType = obj.getType(tmpKey{end});
                % target key
                tmpTargetKey = tmpKey{end};
                % if object leaving efflux particle is pointing to influx node
                if strcmp(tmpType,[myDS.inFluxNode myDS.NODE])
                    % crawl to center
                    tmpKey = obj.linkCrawl(tmpKey{end},2);
                    tmpKey = tmpKey{end};
                    % get center type
                    tmpType = obj.getType(tmpKey);
                    % default to non-matter - assume link crawl is done
                    tmpLinkType = 'nonmatter';
                    % target key
                    tmpTargetKey = tmpKey;
                    % if type is NDML
                    if strcmp(tmpType,[myDS.NDML myDS.NU myDS.NODE])
                        % get the efflux node of NDML
                        tmpKey = obj.getEffluxNode(tmpKey);
                        % get the pointer leaving NDML
                        tmpKey = obj.search({'source_key','type'},{tmpKey,myDS.LINK});
                        tmpKey = tmpKey.key;
                        % crawl to its target - assume isBound
                        tmpKey = obj.linkCrawl(tmpKey,3);
                        tmpKey = tmpKey{end};
                        % 
                        tmpType = obj.getType(tmpKey);
                        tmpLinkType = 'matter';
                        % target key
                        tmpTargetKey = tmpKey;
                    end
                end
                pointerStruct(p).sourceKey = nuKey;
                pointerStruct(p).Type = tmpType;
                pointerStruct(p).linkType = tmpLinkType;
                pointerStruct(p).targetKey = tmpTargetKey;
            end
            pointerStruct = orderfields(pointerStruct,[5 1 2 4 3]);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put object to table
        % ----------------------------
        % purpose: render an en particle to database
        % data is the data for the -edge-node (en)
        % type is the type string for the en
        % target is the link to the target - in the case of direct interaction
        % source is the link to the source - in the case of direction interation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnKey] = put(obj,data,type,target_key,source_key)
            % if there is no source link
            if nargin == 4
                source_key = [];
            end
            % generate time-key
            [K,T] = obj.generateKeyTime();
            % perform insert
            %mksqlite('INSERT INTO flinks VALUES (?,?,?,?,?,?)', {K,T,data,type,target_key,source_key});
            mksqlite('INSERT INTO flinks VALUES (?,?,?,?,?,?,?)', {'',K,T,data,type,target_key,source_key});
            % return key
            returnKey = K;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get object/link from table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [q] = get(obj,key,deepLevel)
            q = mksqlite(['SELECT * FROM flinks WHERE key=''' key '''']);
            if deepLevel >= 1
                if contains(q.type,'_nu_node')
                    deltaKey = obj.getDeltaNode(key);
                    piKey = obj.getPiNode(key);
                    lambdaKey = obj.getLambdaNode(key);
                    efKey = obj.getEffluxNode(key);
                    inKey = obj.getInfluxNode(key);
                    q.deltaNode = obj.get(deltaKey,0);
                    q.piNode = obj.get(piKey,0);
                    q.lambdaNode = obj.get(lambdaKey,0);
                    q.efNode = obj.get(efKey,0);
                    q.inNode = obj.get(inKey,0);
                end
            end
            if deepLevel >=2
                [keys] = obj.getFlagKeys(key);
                for k = 1:numel(keys)
                    q.piNode.flags.(keys{k}) = obj.getFlag(key,keys{k});
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put node to table
        %%%%%
        % data :- data to put into the node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnKey] = putNode(obj,data,nodeType,linkKey)
            % append _ to the node type if needed
            if ~isempty(nodeType)
                if ~strcmp(nodeType(1),'_')
                    nodeType = ['_' nodeType];
                end
            end
            % attach the node type to the string in myDS.NODE to say this is a node
            nodeType = [nodeType myDS.NODE];
            % if the node has a direct link - then render in linkKey
            returnKey = obj.put(data,nodeType,linkKey);
            % render to adj table
            if ~isempty(linkKey) && ~isempty(returnKey)
                obj.addAdjTable(returnKey,'',linkKey,'direct',1);
            end
        end
        
        
        
        
        
        function [] = addAdjTable(obj,sK,lK,tK,renderType,N)
            sK = obj.getRowid(sK);
            tK = obj.getRowid(tK);
            if strcmp(renderType,'links')
                lK = obj.getRowid(lK);
                obj.adjTables{N}(sK,lK) = 1;
                obj.adjTables{N}(lK,tK) = 1;
            else
                obj.adjTables{N}(sK,tK) = 1;
            end
        end
        
        function [] = remAdjTable(obj,sK,tK,N)
            sK = obj.getRowid(sK);
            tK = obj.getRowid(tK);
            obj.adjTables{N}(sK,tK) = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put link to table - quick format syntax - non-matter only
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [linkKey] = link(obj,sourceKey,targetKey,modType)
            isForced = true;
            isDirect = false;
            isMatter = false;
            verbose = true;
            isBidirectional = false;
            bindType = {'ef','in'};
            [linkKey,worked] = obj.putLink(isBidirectional,isForced,isDirect,isMatter,bindType,sourceKey,targetKey,verbose);
            if ~isempty(modType)
                if ~strcmp(modType(1),'_')
                    modType = ['_' modType];
                end
                obj.modifyType(true,linkKey,modType);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put link to table - long format syntax
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [linkKey,worked] = putLink(obj,isBidirectional,isForced,isDirect,isMatter,bindType,sourceKey,targetKey,verbose)
            if verbose
                myDS.generateLinkingCommand(isForced,isDirect,isMatter,sourceKey,targetKey,verbose);
                if isBidirectional
                    myDS.generateLinkingCommand(isForced,isDirect,isMatter,targetKey,sourceKey,verbose);
                end
            end
            
            
            %isBound = true;
            
            if isDirect
                linkKey = sourceKey;
                [worked] = updateDirectLink(obj,isForced,sourceKey,targetKey);
            else
                if isMatter
                    [linkKey,worked] = obj.upsertIndirectMatterLink(isForced,bindType,sourceKey,targetKey);
                else
                    [linkKey,worked] = obj.upsertIndirectNonMatterLink(isForced,bindType,sourceKey,targetKey);
                end
            end
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % link types - direct, non-direct-non-matter, non-direct-matter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [worked] = updateDirectLink(obj,isForced,key,targetKey)
            if obj.keyExists(key)
                if isForced
                    sql = ['UPDATE flinks SET target_key = ''' targetKey ''' WHERE key=''' key ''''];
                    q = mksqlite(sql);
                    worked = true;
                else
                    if ~obj.isDirectLinked(key)
                        sql = ['UPDATE flinks SET target_key = ''' targetKey ''' WHERE key=''' key ''''];
                        q = mksqlite(sql);
                        worked = true;
                    else
                        fprintf('Direct link UPDATE request failed due to existing link.\n');
                        worked = false;
                    end 
                end
            else
                fprintf(['Error: source key does not exisits in this database\n']);
                worked = false;
            end
            % render to changes adjTable
            if worked
                fprintf(['DONT USE YET.\n']);
            end
        end
        
        function [linkKey,worked] = upsertIndirectNonMatterLink(obj,isForced,bindType,sourceKey,targetKey)
            if ~obj.keyExists(sourceKey)
                fprintf(['Warning: Issuing link from source outside the scope of current data store.\n'])
            end
            if ~obj.keyExists(targetKey)
                fprintf(['Warning: Issuing link to target outside the scope of current data store.\n'])
            end
            % non-direct non-matter link type
            linkType = [myDS.LINK];
            
            
            originalSourceKey = sourceKey;
            originalTargetKey = targetKey;
            
            % if isBound - then get the ports/receptors to bind
            if ~isempty(bindType)
                sourceKey = obj.getXNode(sourceKey,bindType{1});
                targetKey = obj.getXNode(targetKey,bindType{2});
            end
            
            
            if isForced
                linkKey = obj.put([],linkType,targetKey,sourceKey);
            else
                if ~(obj.isNonDirectNonMatterLinked(sourceKey,targetKey))
                    linkKey = obj.put([],linkType,targetKey,sourceKey);
                else
                    fprintf('Non-Direct Non-Matter link UPDATE request failed due to existing link.\n');
                    linkKey = [];
                end
            end
            worked = ~isempty(linkKey);
            % add to link table
            if worked
                obj.addAdjTable(originalSourceKey,linkKey,originalTargetKey,'links',2);
            end
        end
        
        function [linkKey,worked] = upsertIndirectMatterLink(obj,isForced,bindType,sourceKey,targetKey)
            if ~obj.keyExists(sourceKey)
                fprintf(['Warning: Issuing link from source outside the scope of current data store.\n'])
            end
            if ~obj.keyExists(targetKey)
                fprintf(['Warning: Issuing link to target outside the scope of current data store.\n'])
            end
            % non-direct non-matter link type
            linkType = [myDS.LINK];
            
            originalSourceKey = sourceKey;
            originalTargetKey = targetKey;
            
            % if isBound - then get the ports/receptors to bind
            if ~isempty(bindType)
                sourceKey = obj.getXNode(sourceKey,bindType{1});
                targetKey = obj.getXNode(targetKey,bindType{2});
            end
            
            
            if isForced
                linkKey = obj.generateNew('','ndml','');
                
                in = obj.getInfluxNode(linkKey);
                out = obj.getEffluxNode(linkKey);
                
                obj.upsertIndirectNonMatterLink(true,'',sourceKey,in);
                obj.upsertIndirectNonMatterLink(true,'',out,targetKey);
            else
                if ~(obj.isNonDirectMatterLinked(sourceKey,targetKey))
                    linkKey = obj.generateNew('','ndml','');
                
                    in = obj.getInfluxNode(linkKey);
                    out = obj.getEffluxNode(linkKey);

                    obj.upsertIndirectNonMatterLink(true,'',sourceKey,in);
                    obj.upsertIndirectNonMatterLink(true,'',out,targetKey);
                else
                    fprintf('Non-Direct Non-Matter link UPDATE request failed due to existing link.\n');
                    linkKey = [];
                end
            end
            worked = ~isempty(linkKey);



            % add to link table
            if worked
                
                obj.addAdjTable(originalSourceKey,linkKey,originalTargetKey,'links',3);
            
            end
        end
        
        
        
        
        
        
        
        function [ret] = isDirectLinked(obj,key)
            ret = ~isempty(obj.queryDirectLink(key));
        end
        
        function [ret] = isNonDirectNonMatterLinked(obj,sourceKey,targetKey)
             ret = ~isempty(obj.queryNonDirectNonMatterLink(sourceKey,targetKey));
        end
        
        function [ret] = isNonDirectMatterLinked(obj,sourceKey,targetKey)
             ret = ~isempty(obj.queryNonDirectMatterLink(sourceKey,targetKey));
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check if key exists
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ret] = keyExists(obj,key)
             sql = ['SELECT key FROM flinks WHERE key = ''' key ''''];
             q = mksqlite(sql);
             ret = ~isempty(q);
        end

        function [exists,key] = typeExists(obj,type)
            key = obj.getTypeKey(type);
            if ~isempty(key)
                exists = true;
            else
                exists = false;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % type key is 1 to 1 if unique flag is on
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [k] = getTypeKey(obj,type)
            type = obj.saneCheck(type);
            q = obj.search({'type'},{[type myDS.NU myDS.NODE]});
            if ~isempty(q)
                k = q.key;
            else
                k = '';
            end
        end

        function [] = putReifiedTypeNode(obj,type,data,varargin)
            what = 1;
        end

        function [props,pathEvidance] = getPotentialProperties(obj,targetKey)
            q = obj.search({'type'},{'_hasa_link'});
            if ~isempty(q)
                q = struct2cell(q);
               
                pathEvidance = obj.dijkstra(q,targetKey,2);
                for e = 1:numel(pathEvidance)
                    tmp = obj.linkCrawlReverse(q{e},1);
                    pathEvidance{e} = [tmp{end} pathEvidance{e}];
                    props{e} = tmp{end};
                end
            else
                pathEvidance = '';
            end
        end

        function [] = proposeProperties(obj,targetKey)
        
        end


        function [path] = getParentProperties(obj,targetKey)
            q = obj.search({'type'},{'_isa_link'});
            q = struct2cell(q);
            for e = 1:numel(q)
                tmp = obj.linkCrawlReverse(q{e},1);
                q{e} = tmp{end};
            end
            path = obj.dijkstra(q,targetKey,2);
        end

        function [path] = getParentTypes(obj,targetKey)
            q = obj.search({'type'},{'_typeof_link'});
            q = struct2cell(q);
            for e = 1:numel(q)
                tmp = obj.linkCrawlReverse(q{e},1);
                q{e} = tmp{end};
            end
            path = obj.dijkstra(q,targetKey,2);
        end





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put a type node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [k] = putTypeNode(obj,modType)
            [exists,k] = obj.typeExists(modType);
            if obj.nonDuplicate && ~exists
                if ~isempty(modType)
                    if ~strcmp(modType(1),'_')
                        modType = ['_' modType];
                    end
                end
                k = obj.generateNew('','tau','',modType);
            end
        end
        
        function [ret] = detailPath(obj,path)
            if ~isempty(path)
                for p = 1:numel(path)
                    type = {};
                    for n = 1:numel(path{p})
                        type{n} = obj.getType(path{p}{n});
                    end
                    ret{p} = cell2table([path{p}' type'],'VariableNames',{'key','type'});
                end
            else
                ret = table();
            end
        end


        function [ret] = detailNodeList(obj,list)
            type = {};
            for n = 1:numel(list)
                type{n} = obj.getType(list{n});
            end
            ret = cell2table([list' type'],'VariableNames',{'key','type'});
        end

        
        function [id] = getRowid(obj,key)
            sql = ['SELECT id FROM flinks WHERE key = ''' key ''''];
            q = mksqlite(sql);
            id = q.id;
        end
        
        function [key] = getKey(obj,id)
            sql = ['SELECT key FROM flinks WHERE id = ' num2str(id)];
            q = mksqlite(sql);
            key = q.key;
        end
        
        function [targetKey] = queryDirectLink(obj,key)
            sql = ['SELECT target_key FROM flinks WHERE key = ''' key ''''];
            q = mksqlite(sql);
            targetKey = q.target_key;
        end
        
        function [linkKey] = queryNonDirectNonMatterLink(obj,sourceKey,targetKey)
            linkKey = obj.search({'source_key','target_key','type'},{sourceKey,targetKey,myDS.LINK});
        end
        
        function [linkKey] = queryNonDirectMatterLink(obj,sourceKey,targetKey)
            [linkKey] = obj.getEffluxLinkList(sourceKey);
            linkKey = linkKey(strcmp({linkKey.targetKey},targetKey));
        end
        
       
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % particle generator - atom too
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [key] = generateNew(obj,data,type,directLink,modType)
            switch type
                case 'ndml'
                    [key] = obj.generateNew([],'tau');
                    obj.modifyType(true,key,myDS.NDML);
                case 'fk'
                    [key] = obj.putNode(data,myDS.FLAGKEY,directLink);
                case 'fv'
                    [key] = obj.putNode(data,myDS.FLAGVALUE,directLink);
                case 'nu'
                    [key] = obj.putNode([],myDS.NU,[]);
                case 'delta'
                    [key] = obj.putNode(data,myDS.DELTA,directLink);
                case 'lambda'
                    [key] = obj.putNode([],myDS.LAMBDA,directLink);
                    [ikey] = obj.putNode([],myDS.inFluxNode,key);
                    [ekey] = obj.putNode([],myDS.efFluxNode,key);
                case 'pi'
                    [key] = obj.putNode([],myDS.PI,directLink);
                case 'tau'
                    [nu_key] = obj.generateNew([],'nu');
                    [delta_key] = obj.generateNew(data,'delta',nu_key);
                    [lambda_key] = obj.generateNew([],'lambda',nu_key);
                    [pi_key] = obj.generateNew([],'pi',nu_key);
                    key = nu_key;
            end
            if nargin == 5
                obj.modifyType(true,key,modType);
            end
        end
        



        function [ret] = modifyType(obj,toCat,key,newType)
            if obj.keyExists(key)
                if toCat
                    newType = [newType obj.getType(key)];
                end
                sql = ['UPDATE flinks SET type = ''' newType ''' WHERE key=''' key ''''];
                q = mksqlite(sql);
                ret = true;
            else
                fprintf(['Error: source key does not exisits in this database.\n']);
                ret = false;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get particle type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [type] = getType(obj,key)
            if obj.keyExists(key)
                sql = ['SELECT type FROM flinks WHERE key=''' key ''''];
                q = mksqlite(sql);
                type = q.type;
            else
                fprintf(['Error: source key does not exisits in this database.\n']);
                type = '';
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put particle data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [k] = putData(obj,data,modType)
            if nargin == 2
                k = obj.generateNew(data,'tau','');
            else
                if ~isempty(modType)
                    if ~strcmp(modType(1),'_')
                        modType = ['_' modType];
                    end
                end
                k = obj.generateNew(data,'tau','',modType);
            end
        end
        

        function [k,worked,oldData] = modData(obj,key,newData)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get particle data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [d] = getData(obj,key)
            delta_key = obj.search({'target_key','type'},{key,[myDS.DELTA myDS.NODE]});
            sql = ['SELECT data FROM flinks WHERE key=''' delta_key.key ''''];
            q = mksqlite(sql);
            d = q.data;
        end
        
        
        function [type] = saneCheck(obj,type)
            if ~isempty(type)
                if ~strcmp(type(1),'_')
                        type = ['_' type];
                end
            end
        end
        
        function [linkKey] = propLink(obj,sourceKey,targetKey,modType)
            isForced = true;
            isDirect = false;
            isMatter = false;
            verbose = true;
            isBidirectional = false;
            bindType = {'nu','in'};
            [linkKey,worked] = obj.putLink(isBidirectional,isForced,isDirect,isMatter,bindType,sourceKey,targetKey,verbose);
            % mod the link = if needed
            if ~isempty(modType)
                [modType] = obj.saneCheck(modType);
                obj.modifyType(true,linkKey,modType);
            end
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flaggin code - mod, add, search are below
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [flagKey] = setFlag(obj,key,flagKey,value)
            piKey = obj.getPiNode(key);
            flagKey = obj.generateNew(flagKey,'fk',piKey);
            valueKey = obj.generateNew(value,'fv',flagKey);
        end
        
        function [value] = getFlag(obj,key,flagKey)
            piKey = obj.getPiNode(key);
            flagKey_key = obj.search({'target_key','type','data'},{piKey,[myDS.FLAGKEY myDS.NODE],flagKey});
            flagValue_key = obj.search({'target_key','type'},{flagKey_key.key,[myDS.FLAGVALUE myDS.NODE]});
            sql = ['SELECT data FROM flinks WHERE key=''' flagValue_key.key ''''];
            q = mksqlite(sql);
            value = q.data;
        end
        
        function [keys] = getFlagKeys(obj,key)
            piKey = obj.getPiNode(key);
            flagKeys = obj.search({'target_key','type'},{piKey,[myDS.FLAGKEY myDS.NODE]});
            for k = 1:numel(flagKeys)
                sql = ['SELECT data FROM flinks WHERE key=''' flagKeys(k).key ''''];
                q = mksqlite(sql);
                keys{k} = q.data;
            end
        end
        
        function [nuKeys] = flagSearchByKey(obj,flagKey)
            res = obj.search({'type','data'},{[myDS.FLAGKEY myDS.NODE],flagKey});
            for k = 1:numel(res)
                [keyChain] = linkCrawl(obj,res(k).key,2);
                nuKeys{k} = keyChain{end};
            end
        end
        
        function [nuKeys] = flagSearchByValue(obj,flagValue)
            res = obj.search({'type','data'},{[myDS.FLAGVALUE myDS.NODE],flagValue});
            for k = 1:numel(res)
                [keyChain] = linkCrawl(obj,res(k).key,3);
                nuKeys{k} = keyChain{end};
            end
        end
        
        
        
        
        
        
        
        
        
        function [keyChain] = linkCrawl(obj,startKey,N)
            keyChain{1} = startKey;
            for e = 1:N
                sql = ['SELECT target_key FROM flinks WHERE key=''' keyChain{e} ''''];
                q = mksqlite(sql);
                keyChain{e+1} = q.target_key;
            end
        end
        
        
        
        
        function [keyChain] = linkCrawlReverse(obj,startKey,N)
            keyChain{1} = startKey;
            for e = 1:N
                sql = ['SELECT source_key FROM flinks WHERE key=''' keyChain{e} ''''];
                q = mksqlite(sql);
                keyChain{e+1} = q.source_key;
            end
        end
        
        
        
        
        function [] = print(obj)
            q = mksqlite( 'SELECT * FROM flinks' );
            q = struct2table(q);
            q
        end
        %function [] = propSearchByKey(obj,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get object/link from table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnTable] = get(obj,keyList)
            for e = 1:numel(keyList)
                IDX(e) = find(obj.ds.key == keyList(e));
            end
            returnTable = obj.ds(IDX,:);
        end
        %}
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % modify object/link from table - ????
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = modify(obj,key,field,value)
            index = obj.getIndex(key);
            if ~isempty(index)
                if strcmp(field,'data')
                    value = {value};
                end
                obj.ds{index,field} = value;
            end
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % basic: init,put(s),get(s),modify(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test if object is a type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [tf] = isType(obj,key,queryType)
            index = obj.getIndex(key);
            tf = strcmp(obj.ds.type(index),queryType);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test if object is a sub-type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [tf] = isSubType(obj,key,queryType)
            index = obj.getIndex(key);
            tf = contains(obj.ds.type(index),queryType);
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get object/link from table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [str] = toString(obj,verbose)
            if nargin == 1
                verbose = false;
            end
            for e = 1:size(obj.ds,1)
                key = obj.ds.key(e);
                str{e} = obj.keyToString(key);
                fprintf(str{e});
                 
                if verbose
                    if obj.isSubType(key,'_queue')
                        obj.queueToString(key);
                    end
                    if obj.isSubType(key,'_context')
                        obj.contextToString(key);
                    end
                end
                
               
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create octet
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [str] = keyToString(obj,key)
            
            index = obj.getIndex(key);
            o1 = num2str(obj.BASE*double(bitget(key,25:32,'uint32'))');
            o2 = num2str(obj.BASE*double(bitget(key,17:24,'uint32'))');
            o3 = num2str(obj.BASE*double(bitget(key,9:16,'uint32'))');
            o4 = num2str(obj.BASE*double(bitget(key,1:8,'uint32'))');
            str = [num2str(key) '\t-->\t' o1 '.' o2 '.' o3 '.' o4 '\t-->\t' obj.ds.type{index} '\n'];
                
        end
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set empty pointer
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = setEmptyPtr(obj)
            obj.emptyPtr = size(obj.ds,1) + 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get pointer to table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ptr] = getPtr(obj)
            ptr = obj.emptyPtr;
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mass import files from list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = importFromTable(obj,pTable)
            for r = 1:size(pTable,1)
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mass import files from list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = importImagesFromList(obj,fileList,type,toQueue)
            for f = 1:numel(fileList)
                fprintf(['Importing fileName:' fileList{f} '.\n']);
                imageKey = obj.putData(fileList{f},type);
                frameKey = obj.putData(eye(3),'refFrame');
                obj.link(imageKey,frameKey,'hasa');
                if toQueue
                    obj.sequence{end+1} = imageKey;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % search for objects by type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnKey] = searchByType(obj,type)
            if ~iscell(type)
                type = {type};
            end
            indexMask = false(size(obj,1),1);
            for e = 1:numel(type)
                indexMask = indexMask | strcmp(obj.ds.type,type{e});
            end
            returnKey = obj.ds.key(indexMask);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % filter objects by type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [keyList] = filterByType(obj,keyList,filterType)
            tf = obj.isType(keyList,filterType);
            keyList = keyList(tf);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % filter objects by type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [keyList] = filterBySubType(obj,keyList,filterType)
            tf = obj.isSubType(keyList,filterType);
            keyList = keyList(tf);
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % search for objects by type and linkkey
        % find objects that have their linkKey as
        % given by queryLinkKey AND have type queryType
        % ----------------------------------
        % find those objects that directly point to an
        % object and are of a type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnKey] = searchByTypelinkKey(obj,queryType,queryLinkKey,linkDirection)
            
            % look for objects that either
            % directly point to(-1) or from(+1)
            if linkDirection < 0
                IDX = find(obj.ds.linkKey == queryLinkKey);
            else
                IDX = find(obj.ds.slinkKey == queryLinkKey);
            end
            
            
            if ~iscell(queryType)
                queryType = {queryType};
            end
            
            indexMask = false(size(IDX,1),1);
            for e = 1:numel(queryType)
                indexMask = indexMask | strcmp(obj.ds.type(IDX),queryType{e});
            end
            
            returnKey = obj.getKey(IDX(indexMask));
        end
        
        
        
        
        
        
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % recursive search
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnKeyList] = recursiveSearch(obj,key)
            % find object
            [index] = getIndex(obj,queryKey);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init sequence by key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = initSequenceByKey(obj,keyList)
            obj.sequence = keyList;
            obj.seqPtr = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init sequence by type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = initSequenceByType(obj,type)
            [returnKey] = searchByType(obj,type);
            obj.initSequenceByKey(returnKey);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init execute sequence for type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [contextKey] = initExecuteSequenceByType(obj,type,executeType)
            [returnKey] = searchByType(obj,type);
            for e = 1:numel(returnKey)
                contextKey(e) = initExecuteContext(obj,returnKey(e),executeType);
            end
            obj.initSequenceByType('executeContext');
        end
        %}
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read N records from the sequence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [data] = readN(obj,N)
            for e = 1:N
                if obj.hasData()
                    data{e} = obj.read();
                else
                    break
                end
            end
        end
        
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of objects in the sequence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [N] = numToRead(obj)
            N = numel(obj.sequence);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the geoframe for the key list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [frameKeyList] = getFrameForKeyList(obj,keyList)
            [index] = obj.getIndex(keyList);
            for e = 1:numel(index)
                frameKeyList(e) = obj.ds.linkKey(index(e));
            end 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call crop - write geoframe and write image 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = imcrop(obj,sourceImage,cropBox,objectType,contextKey)
            tmpI = imcrop(sourceImage,cropBox);
            
        end
        
        function [] = imwrite(obj,image,objectType,contextKey,frameKey)
            
            % cat _image for base class of image
            objectType = [objectType '_image'];
            
            
            
            % get the context type
            context_type = obj.ds.data{obj.getIndex(contextKey)};
            % get the key for the context
            contextState = obj.getExecuteContextKey(contextKey,'context_state');
            
            
            % generate random name
            [fName,~] = obj.generateKeyTime();
            % convert number to string
            fName = num2str(fName);
            % hash value
            hashValue = [fName(1:2) filesep];
            % construct full object type
            fullObjectType = [context_type '_' objectType];
            
            
            
            
            
            % sort for remote last
            for f = 1:numel(obj.base_outPath)
                sortValue(e) = uint8(isIRODS(obj.base_outPath{e}));
            end
            [sortValue,sidx] = sort(sortValue);
            obj.base_outPath = obj.base_outPath(sidx);
            
            
            % loop and write
            for f = 1:numel(obj.base_outPath)
                localName = '';
                
                % construct outpath and file name
                outPath{f} = [obj.base_outPath{f} context_type filesep objectType filesep hashValue filesep];
                fileName{f} = [outPath{f} fName '.tif'];

                
                if isIRODS(outPath{f})
                    CMD = ['imkdir -p ' outPath];
                else
                    CMD = ['mkdir -p ' outPath];
                end
                system(CMD);

                
                if isIRODS(outPath{f})
                    
                    obj.iimwrite(image,fileName);
                else
                    localName = fileName;
                    imwrite(image,fileName);
                end
            end
            
            
           
           
         
            
            
            obj.put(fileName,fullObjectType,frameKey);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % irods interface
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get ticket from context
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ticket] = getTicket(obj,imageKey)
            ticket = '';
            if isdeployed()
                if isIRODS(imageKey)
                    fidx = strfind(imageKey,filesep);
                    userName = imageKey((fidx(3)+1):(fidx(4)-1));
                    ticket = getTicket(userName,'read','./private.key','./');
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push to irods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = ipush(obj,localFileList,remoteFileList,ticket)
            
            if ~iscell(localFileList)
                localFileList{1} = localFileList;
            end
            if ~iscell(remoteFileList)
                remoteFileList{1} = remoteFileList;
            end
            
            if ~isempty(location)
                fprintf(['**************************************************************************\n'])
                fprintf(['Start push to irods:' num2str(numel(localFileList)) ' files.\n']);tic
                fprintf(['**************************************************************************\n'])
                for e = 1:numel(localFileList)
                    fprintf(['*************************************\n']);
                    fprintf(['Start push to irods:' num2str(e) ':' num2str(numel(localFileList)) ' files.\n']);tic
                    fprintf(['*************************************\n']);
                    [~,tn,te] = fileparts(localFileList{e});
                    targetFile = [pth filesep tn te];
                    targetHTTP = [targetFile '#' ticket '#'];
                    targetHTTP = xform2URL({targetHTTP});
                    targetHTTP = targetHTTP{1};
                    cmd = ['iput -f -V -t ' ticket ' "' localFileList{e} '" "' targetFile '"'];
                    %cmd = ['curl -X POST ' targetHTTP ' -F uploadFile=@' fileList{e}];
                    fprintf(['Push Command is:' cmd '\n']);
                    [r,o] = system(cmd,'-echo');
                    fprintf(['\npushing to file to irods:' tn '\n']);
                    fprintf(['*************************************\n']);
                    fprintf(['End push to irods:' num2str(e) ':' num2str(numel(localFileList)) ' files.\n']);tic
                    fprintf(['*************************************\n']);
                end
                fprintf(['**************************************************************************\n'])
                fprintf(['End push to irods:' num2str(toc) ' seconds.\n']);
                fprintf(['**************************************************************************\n'])
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % push to irods
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %{
        % raw iget
        function [] = iget(obj,file)
        
        end
        
        % raw iput
        function [] = iput(obj,file)
        
        end
        %}
        % irods imwrite 
        function [] = iimwrite(obj,image,remoteName)
            % create local file name
            localTmpName = tempname;
            localTmpName = [localTmpName(6:end) '.tif'];
            % write local file
            imwrite(image,localTmpName);
            %  
            [remotePath,nm,ext] = fileparts(remoteName);
            % 
            pushToiRods(rPath,fileList);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
       
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % methods for structures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SEQUENCE START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SEQUENCE END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONTEXT START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init context type object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [contextKey] = initContext(obj,contextType,initFKV,initEFL)
            contextType = [contextType myDS.CONTEXT];
            % create execute context
            contextKey = obj.putNode('',contextType,'');


            for e = 1:numel(initFKV)
                % init the request time
                fieldRecord(e) = obj.initField(initFKV(e).value,initFKV(e).key,contextKey);
            end

            for e = 1:numel(initEFL)
                linkRecord(e) = obj.putLink(contextKey,initEFL(e).type,initEFL(e).key);
            end
            
        end

        function [contextKey] = initExecuteContext(obj,fileKey)
        
            contextType = '_execute';
            initFKV = [];
            initFKV(1).key ='init_time';
            initFKV(1).value = datenum(datetime);
            initFKV(2).key ='context_state';
            initFKV(2).value = 0;
            initFKV(3).key ='start_execute_time';
            initFKV(3).value = 0;

            initFKV(4).key ='end_execute_time';
            initFKV(4).value = 0;



            initEFL(1).type = 'file';
            initEFL(1).key = fileKey;

            [contextKey] = obj.initContext(contextType,initFKV,initEFL);
        end

        function [] = initAlgoContext(obj,func)
            
            varSpace = functions(func);
            S = varSpace.workspace{1};
           
            fStr = func2str(func);
            s1 = strfind(fStr,'(');
            s2 = strfind(fStr,')');
            matlabFunction = fStr(s2(1)+1:s1(2)-1);


            flds = fields(varSpace.workspace{1});
            for f = 1:numel(flds)
                initFKV(f).key = flds{f};
                initFKV(f).value = S.(flds{f});
            end
        
            

            here = 1;
        end

        function [] = contextToString(obj,contextKey)
            sep = ['\t\t\t----------------------------------------------------------------------------------\n'];
            keyString = obj.keyToString(contextKey);
            keyString((end-1):end) = [];
            str{1} = ['\t\t\t|MN:' keyString  '\t | \n'];
            
            fprintf(sep);
            fprintf(str{1});
            fprintf(sep);

            fields = obj.getFieldKeys(contextKey);
            for e = 1:numel(fields)
                keyString = obj.keyToString(fields(e));
                keyString((end-1):end) = [];
                fldstr{e} = ['\t\t\t|FLD:' keyString '\t | \n'];
                fprintf(fldstr{1});
            end
            fprintf(sep);
            links = obj.getFlux(contextKey,8);
            for e = 1:numel(links)
                keyString = obj.keyToString(links(e));
                keyString((end-1):end) = [];
                lkstr{e} = ['\t\t\t|LNK:' keyString '\t | \n'];
                fprintf(lkstr{1});
            end
            fprintf(sep);


        end

        function [] = executeFunctionOnKeySequence(obj,keySet,func)
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mod context fields
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [contextInfo,fileKey,algoKey] = getContext(obj,contextKey)

            % init the algo key to nothing
            algoKey = '';
            %{
            % find all index of objects under the context key
            IDX = find(obj.ds.linkKey == contextKey);
            % get keys of index
            KEY = obj.getKey(IDX);

            % find the file link for the context
            fidx = find(strcmp(obj.ds.type(IDX),'file_context_link'));
            % get the file key
            fileKey = obj.ds.data{IDX(fidx)};
            % get the file index
            fileIndex = obj.getIndex(fileKey);
            % get the frame key
            frameKey = obj.ds.linkKey(fileIndex);
            contextKeys = KEY;
            %}
            % gather context links
            linkSet = obj.getFlux(contextKey,[1+8]);

            % stack keys
            contextKeys = [contextKeys;fileKey;frameKey];
            contextIDX = obj.getIndex(contextKeys);
            contextInfo = obj.ds(contextIDX,:);

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mod context fields
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = modContext(obj,modType,data,contextKey)
            switch modType
                case 'state'
                    returnKey = obj.getExecuteContextKey(contextKey,'context_state');
                case 'startTime'
                    returnKey = obj.getExecuteContextKey(contextKey,'start_execute_time_stamp');
                case 'endTime'
                    returnKey = obj.getExecuteContextKey(contextKey,'end_execute_time_stamp');
                case 'filexFertime'
            end
            obj.modField(returnKey,data);
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONTEXT START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIELD START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init field type object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fieldKey] = initField(obj,data,fieldType,parentKey)
            fieldType = [fieldType '_field'];
            % create execute context
            fieldKey = obj.putNode(data,fieldType,parentKey,true);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mod field data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fieldKey] = modField(obj,fieldKey,data)
            index = obj.getIndex(fieldKey);
            obj.ds.data{index} = {data};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get field keys for parent key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [inFlux] = getFieldKeys(obj,parentKey,subType)
            inFlux = obj.getFlux(parentKey,1);
            tf = obj.isSubType(inFlux,'_field');
            inFlux = inFlux(tf);
            if nargin == 3
                tf = obj.isSubType(inFlux,subType);
                inFlux = inFlux(tf);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIELD START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % QUEUE START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init queue type object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [queueKey] = initQueue(obj,queueType)
            queueType = [queueType myDS.QUEUE];
            queueKey = obj.putNode('',queueType,'');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init queue node type object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [newQueueElement] = initQueNode(obj,queueData,queueType)
            queueType = [queueType myDS.QUEUEELE];
            newQueueElement = obj.putNode(queueData,queueType,'');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init queue type object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = readQueue(obj,queueKey,queuePtrKey)
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % is queue empty
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [tf] = isQueueEmpty(obj,queueKey)
            [targetKeyList] = obj.getEffluxLinkList(queueKey);
            tf = isempty(targetKeyList);
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the queue tail key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [tailKey] = getQueueTail(obj,queueKey)
            tailKey = [];
            if ~obj.isQueueEmpty(queueKey)
                EffluxLinkKeyList = obj.getEffluxLinkList(queueKey);
                if ~isempty(EffluxLinkKeyList)
                    tailKey = obj.filterByType(EffluxLinkKeyList,'queueTail_link');
                end
                tailKey = obj.traverseLink(tailKey,'forward');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the queue head key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [headKey] = getQueueHead(obj,queueKey)
            headKey = [];
            if ~obj.isQueueEmpty(queueKey)
                EffluxLinkKeyList = obj.getEffluxLinkList(queueKey);
                if ~isempty(EffluxLinkKeyList)
                    headKey = obj.filterByType(EffluxLinkKeyList,'queueHead_link');
                end
                headKey = obj.traverseLink(headKey,'forward');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assign queue head
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = assignQueueHead(obj,queueKey,headKey)
           obj.putLink(queueKey,'queueHead',headKey);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assign queue tail
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = assignQueueTail(obj,queueKey,tailKey)
            obj.putLink(queueKey,'queueTail',tailKey);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push data to queue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = pushToDataQueue(obj,queueKey,queueData,queueType)
            
            
            % create new queue element
            newQueueElement = obj.initQueNode(queueData,queueType);
            
            % if queue is empty
            if obj.isQueueEmpty(queueKey)
                obj.assignQueueHead(queueKey,newQueueElement);
%{
            else
                EffluxLinkKeyList = obj.getEffluxLinkList(queueKey);
                HeadLinkKeyList = obj.filterByType(EffluxLinkKeyList,'queueHead_link');
                obj.modify(HeadLinkKeyList,'queueHead_link',newQueueElement);
%}
            end
            
            % if queue is empty - assign tail else mod tail to new element
            EffluxLinkKeyList = obj.getEffluxLinkList(queueKey);
            TailLinkKeyList = obj.filterByType(EffluxLinkKeyList,'queueTail_link');
           
            if isempty(TailLinkKeyList)
                obj.assignQueueTail(queueKey,newQueueElement);
            else
                % add on to the tail
                tailElement = obj.traverseLink(TailLinkKeyList,'forward');
                obj.putLink(tailElement,'next',newQueueElement);
                % change thte tail
                obj.modify(TailLinkKeyList,'linkKey',newQueueElement);
            
            end


            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % render queue to string
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = queueToString(obj,queueKey)
            try
                sep = ['----------------------------------------------------------------------------------\n'];
                tailKey = obj.getQueueTail(queueKey);
                headKey = obj.getQueueHead(queueKey);

                if isempty(tailKey)
                    tailString = ['\t\t\t|TL: NaN' '\t-->\t' 'tailPointer \t|\n'];
                else
                    keyString = obj.keyToString(tailKey);
                    keyString((end-1):end) = [];
                    tailString =['\t\t\t|TL:' keyString '\t | \n'];
                end

                if isempty(headKey)
                    headString = ['\t\t\t|HD: NaN' '\t-->\t' 'headPointer \t|\n'];
                else
                    keyString = obj.keyToString(headKey);
                    keyString((end-1):end) = [];
                    headString =  ['\t\t\t|HD:' keyString '\t | \n'];
                end


                if ~isempty(headKey)
                    currentKey = headKey;
                    str ={};
                    cnt = 1;
                    while currentKey ~= tailKey
                        keyString = obj.keyToString(currentKey);
                        keyString((end-1):end) = [];
                        str{cnt} = ['\t\t\t|ND:' keyString '\t | \n'];
                        cnt = cnt + 1;
                        currentKey = obj.next(currentKey);
                    end
                    keyString = obj.keyToString(currentKey);
                    keyString((end-1):end) = [];
                    str{cnt} = ['\t\t\t|ND:' keyString '\t | \n'];
                end


                fprintf(['\t\t\t' sep]);
                fprintf(headString);
                fprintf(tailString);
                fprintf(['\t\t\t' sep]);
                for e = 1:numel(str)
                    fprintf(str{e});
                end
                fprintf(['\t\t\t' sep]);
            catch ME
                here =1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get next from object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [nextKey] = next(obj,queueEle)
            nextKey = [];
            if ~obj.isQueueEmpty(queueEle)
                EffluxLinkKeyList = obj.getEffluxLinkList(queueEle);
                if ~isempty(EffluxLinkKeyList)
                    nextKey = obj.filterByType(EffluxLinkKeyList,'next_link');
                end
                nextKey = obj.traverseLink(nextKey,'forward');
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % QUEUE END
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods for structures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods for datastore
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        function [data,info] = read(obj)
            %%%%%%%%%%%%%%%%%%%%
            % get the key for the "next" object to be read
            %nextKey = obj.sequence(obj.seqPtr);
            %%%%%%%%%%%%%%%%%%%%
            % init info for object to be read
            [info.executeTable,fileKey] = obj.getExecuteContext(nextKey);
            info.fileKey = fileKey;
            info.executeKey = nextKey;
            % get the file index
            fileIndex = obj.getIndex(fileKey);
            info.frameKey = obj.ds.linkKey(fileIndex);
            % get the data - file anem
            data = obj.ds.data{fileIndex};
            % get the ticket if needed
            info.ticket = obj.getTicket(data);
            % increment the pointer
            obj.seqPtr = obj.seqPtr + 1;
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read record from the sequence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [data,key] = read(obj)
            key = obj.sequence{obj.seqPtr};
            data = obj.getData(key);
            if contains(obj.getType(key),'fileName')
                if obj.toReadThrough
                    % change to general reader later
                    data = imread(data);
                end
            end
            obj.seqPtr = obj.seqPtr + 1;
        end
        
        
        function [tf] = hasdata(obj)
            tf = numel(obj.sequence) >= obj.seqPtr;
        end

        function [] = reset(obj)
            obj.seqPtr = 1;
        end

        function [percent] = progress(obj)
            percent = (obj.seqPtr - 1)/numel(obj.sequence);
        end

        function [data] = preview(obj)
            N = 8;
            [data] = obj.readN(obj);
        end

        function [data] = readAll(obj)
            N = obj.numToRead();
            data = obj.readN();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods for datastore
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods(s) for parallel datastore
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [sub] = partition(obj,n,index)
            % if there is data present
            if obj.hasdata()
                % find the proper index to the sequence list
                subIndex = find((rem(1:numel(obj.sequence),n)+1) == index);
                % get the keys for the subset
                subKey = obj.sequence(subIndex);
                % get the associated keys to package
                totalKeyList = [];
                for k = 1:numel(subKey)
                    [contextInfo] = obj.getExecuteContext(subKey(k));
                    totalKeyList = [totalKeyList ;contextInfo.key];
                end
                % get the index of the keys to package
                totalIndex = obj.getIndex(totalKeyList);
                % package and return
                sub = myDS(obj.base_outPath,obj.ds(totalIndex,:),subKey);
            else
                sub = [];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods(s) for parallel datastore
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % QUEUE END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    methods (Access=protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % max partitions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [n] = maxpartitions(obj)
            n = floor(numel(obj.sequence)/obj.filesPerPartitions);
        end
    end

    methods (Access=private)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % basic create and cast key <--> index
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate key and time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [K,T] = generateKeyTime(obj)
            % generate random index key
            K = randi((2^32)-1,1,'uint32');
            %K = typecast([uint32(0) K],'double');
            K = num2str(K);
            %K = int64(K);
            T = datenum(datetime);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transform keyList to indexList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [indexList] = getIndex(obj,keyList)
            % init to no index found
            indexList = NaN*zeros(size(keyList));
            for e = 1:numel(keyList)
                tmpIndex = find(obj.ds.key == keyList(e));
                if ~isempty(tmpIndex)
                    indexList(e) = tmpIndex;
                end
            end
        end
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transform indexList to keyList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [keyList] = getKey(obj,indexList)
            % init to no index found
            if isempty(indexList)
                keyList = NaN;
            else
                keyList = obj.ds.key(indexList);
            end
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % basic create and cast key <--> index
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % used to execute linking commands - NEEDED
        function [] = renderLink(isForced,isDirect,isMatter,sourceKey,targetKey,verbose)
            
            
            if isDirect
                if isForced
                    
                else
                    
                end
            else
                
            end
            
            
            cmd = '';
            switch isForced
                case 0
                    cmd = [cmd 'non-forced'];
                case 1
                    cmd = [cmd 'forced'];
                otherwise
                    cmd = [cmd 'non-forced'];
            end
            
            
            switch isDirect
                case 0
                    cmd = [cmd ' non-direct' ];
                    switch isMatter
                        case 0
                            cmd = [cmd ' non-matter' ];
                        case 1
                            cmd = [cmd ' matter' ];
                        otherwise
                    end
                case 1
                    cmd = [cmd ' direct' ];
                    cmd = [cmd ' non-matter'];
                otherwise
                    cmd = [cmd ' direct' ];
                    cmd = [cmd ' non-matter'];
            end
            
            cmd = [cmd ' ' sourceKey ':-->' targetKey '\n'];
            fprintf(cmd);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate lambda/link node pairs for data/delta node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % render inFlux object for key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = renderInfluxNode(obj,key)
            obj.putNode('',myDS.efFluxNode,key);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % render inFlux object for key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = renderEffluxNode(obj,key)
            obj.putNode('',myDS.inFluxNode,key);
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get flux by type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        function [linkKeys] = getFlux(obj,nodeKey,fluxNumber)
            linkKeys = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the direct influx links
            if bitget(fluxNumber,1)
                fidx = find(obj.ds.linkKey == nodeKey);
                fidxKey = obj.getKey(fidx);
                subIDX = (obj.isSubType(fidxKey,obj.NODE));
                linkKeys = [linkKeys;fidxKey(subIDX)];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the indirect influx links
            if bitget(fluxNumber,2)
                fidx = find(obj.ds.linkKey == nodeKey);
                fidxKey = obj.getKey(fidx);
                subIDX = (obj.isSubType(fidxKey,obj.LINK));
                linkKeys = [linkKeys;fidxKey(subIDX)];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the direct (one) efflux link(s)
            if bitget(fluxNumber,3)
                index = obj.getIndex(nodeKey);
                fidxKey = obj.ds.linkKey(index);
                linkKeys = [linkKeys;fidxKey(subIDX)];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the indirect (one) efflux link(s)
            if bitget(fluxNumber,4)
                efFluxKey = obj.getEffluxNode(nodeKey);
                fidx = find(obj.ds.slinkKey == efFluxKey);
                fidxKey = obj.getKey(fidx);
                subIDX = (obj.isSubType(fidxKey,obj.LINK));
                linkKeys = [linkKeys;fidxKey(subIDX)];
            end     
        end
        %}
        %{
        function [targetKey] = traverseLink(obj,linkKey,direction)
            if strcmp(direction,'forward')
                if obj.isSubType(linkKey,myDS.NODE)
                    index = obj.getIndex(linkKey);
                    targetKey = obj.ds.linkKey(index);
                elseif obj.isSubType(linkKey,myDS.LINK)
                    index = obj.getIndex(linkKey);
                    fluxNode = obj.ds.linkKey(index);
                    index = obj.getIndex(fluxNode);
                    targetKey = obj.ds.linkKey(index);
                end
            elseif strcmp(direction,'reverse')
                if obj.isSubType(linkKey,myDS.NODE)
                    targetKey = linkKey;
                elseif obj.isSubType(linkKey,myDS.LINK)
                    index = obj.getIndex(linkKey);
                    fluxNode = obj.ds.slinkKey(index);
                    index = obj.getIndex(fluxNode);
                    targetKey = obj.ds.linkKey(index);
                end
            end
        end
        %}
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get inFlux object for node keyList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [inFluxKey] = getInfluxNode(obj,keyList)
            for e = 1:numel(keyList)
                [inFluxKey(e)] = obj.searchByTypelinkKey([myDS.inFluxNode myDS.NODE],keyList);
            end
        end
        %}
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get inFlux object for node keyList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [efFluxKey] = getEffluxNode(obj,keyList)
            for e = 1:numel(keyList)
                [efFluxKey(e)] = obj.searchByTypelinkKey([myDS.efFluxNode  myDS.NODE],keyList);
            end
        end
        %}
        
        
        
        
        
            
        
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods for external disk(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scan irods at dataPath for files with ext
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [FileList] = gdig_irods(dataPath,FileExt)
            %dataPath = ['/iplant/home/' user '/' plantType 'Data/' tissueType 'Data%'];
            if strcmp(dataPath(end),filesep)
                dataPath(end) = [];
            end
            if ~strcmp(dataPath(end),'%')
                dataPath = [dataPath '%'];
            end
            CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
            [o,r] = system(CMD);
            [r] = parseRecords(r);
            for e = 1:numel(r)
                FileList{e} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
            end
            if numel(FileExt) ~= 0
                rm = [];
                for e = 1:numel(FileList)
                    [p,n,ext] = fileparts(FileList{e});
                    
                    if ~any(strcmp(lower(ext(2:end)),FileExt))
                        rm = [rm e];
                    end
                end
                FileList(rm) = [];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scan local for file types
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [FileList] = gdig(FilePath,FileExt)
            FileList = gdig(FilePath,{},FileExt,1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % structured local dig
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [FileList] = sdig(FilePath,FileExt)
            FileList = sdig(FilePath,{},FileExt,1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % image dig at irods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [FileList] = igdig_irods(dataPath)
            FileList = myDS.gdig_irods(dataPath,myDS.IMAGE_EXT_LIST);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % meta dig at irods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [FileList] = mgdig_irods(dataPath)
            FileList = myDS.gdig_irods(dataPath,myDS.META_EXT_LIST);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % nlp on fileList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [pathStats,fileStats,totalPathStats,totalFileStats] = keywordParseFileList(FileList,keywordList)
            if nargin == 1
                keywordList = myDS.KEYWORD;
            end
            
            
            for e = 1:numel(FileList)
                [pth,nm] = fileparts(FileList{e});
                for k = 1:numel(myDS.KEYWORD)
                    P_RES(e,k) = uint8(contains(pth,keywordList{k}));
                    F_RES(e,k) = uint8(contains(nm,keywordList{k}));
                end
            end
            pathStats = array2table(P_RES);
            fileStats = array2table(F_RES);
            pathStats.Properties.VariableNames = keywordList;
            fileStats.Properties.VariableNames = keywordList;
            
            
            columNames = pathStats.Properties.VariableNames;
            for e = 1:numel(columNames)
                totalPathStats.(columNames{e}) = sum(pathStats.(columNames{e}));
            end
            
            columNames = fileStats.Properties.VariableNames;
            for e = 1:numel(columNames)
                totalFileStats.(columNames{e}) = sum(fileStats.(columNames{e}));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % nlp on fileList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [bag] = gatherWordBag(FileList)
            bag = {};
            cnt = 1;
            for e = 1:numel(FileList)
                [pth,nm] = fileparts(FileList{e});
                pth = [pth filesep];
                sidx = strfind(pth,filesep);
                for s = 1:(numel(sidx)-1)
                    newWord = pth(sidx(s)+1:sidx(s+1)-1);
                    bag{cnt} = newWord;
                    fprintf(['added:' newWord '\n']);
                    cnt = cnt + 1;
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods for external disk(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods for external disk(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = generateLinkingCommand(isForced,isDirect,isMatter,sourceKey,targetKey,verbose)
            cmd = '';
            
            switch isForced
                case 0
                    cmd = [cmd 'non-forced'];
                case 1
                    cmd = [cmd 'forced'];
                otherwise
                    cmd = [cmd 'non-forced'];
            end
            
            
            switch isDirect
                case 0
                    cmd = [cmd ' non-direct' ];
                    switch isMatter
                        case 0
                            cmd = [cmd ' non-matter' ];
                        case 1
                            cmd = [cmd ' matter' ];
                        otherwise
                    end
                case 1
                    cmd = [cmd ' direct' ];
                    cmd = [cmd ' non-matter'];
                otherwise
                    cmd = [cmd ' direct' ];
                    cmd = [cmd ' non-matter'];
            end
            
            cmd = [cmd ' ' sourceKey ':-->' targetKey '\n'];
            fprintf(cmd);
        end
    end
end

%{


    FilePath = '/mnt/spaldingdata/nate/mirror_images/rue/';
    FileList = {};
    FileExt = {'tif','TIF'};

    dataPath = '/iplant/home/turnersarahd/Arabidopsis_scans/';
    FileList = myDS.gdig_irods(dataPath,{'tiff'});



    dataPath = '/iplant/home/turnersarahd/';
    dataPath = '/iplant/home/kmichel/';
    dataPath = '/iplant/home/hirsc213/maizeData/seedlingData/';


    FileList = myDS.igdig_irods(dataPath);
    [pathStats,fileStats,totalPathStats,totalFileStats] = myDS.keywordParseFileList(FileList);
    [bag] = myDS.gatherWordBag(FileList);


    dataPath = '/iplant/home/kmichel/analyses/';
    metaFileList = myDS.mgdig_irods(dataPath);
    

    dataStore = myDS();
    dataStore.generateNew('tau')
    tID = dataStore.searchByType('_ta')
    [lambdaKey] = dataStore.getTauLinkComplex(tID);
   
%}


%{



        %{
        function [] = putLinkType(obj,linkType)
            % generate time-key
            [K,T] = obj.generateKeyTime();
            mksqlite('INSERT INTO tlinks VALUES (?,?,?,?,?,?,?)', {'',K,T,linkType,parent});
        end
        
        function [ret] = linkTypeExists(obj,linkType)
              sql = ['SELECT key FROM tlinks WHERE  = linkType ''' linkType ''''];
        end
        
        function [] = upsertlinkType(obj,linkType)
           
        end
        
        function [key] = getLinkTypeKey(obj,linkType)
            sql = ['SELECT key FROM tlinks WHERE  = linkType ''' linkType ''''];
            q = mksqlite(sql);
            if ~isempty(q)
                key = q.key;
            end
        end
        
        
        function [parent,parentKey] = parselinkTypeParent(obj,nonMatterLinkType)
            fidx = strfind(nonMatterLinkType,'_');
            if numel(fidx) > 1
                parent = nonMatterLinkType(fidx(2):end);
                parentKey = obj.getLinkTypeKey(linkType);
            else
                parent = '';
                parentKey = '';
            end
        end
        %}
        
%}