classdef myDS < handle & matlab.io.Datastore & matlab.io.datastore.Partitionable
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % datastore: holds nodes and edges
    % this is starting to rebuild the datastore that I need
    % the basis particle is a node/edge object - enode or nedge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % each basic particle has potential for
    %      1: source, target, and data
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
    % the eLambda handles the efflux attachments to other Tau particles
    % the iLambda handles the infllux attachments to other Tau particles
    %%%%%%%
    % delta: this particle will house the data. as simple as that. for now
    %%%%%%%
    % pi: this particle handles the properties for the data/associated node
    
    
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % main data store objects
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data store
        ds;
        % current read location
        readPtr;
        % current insert location
        emptyPtr;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sequence path through the tree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sequence
        sequence;
        % sequence pointer
        seqPtr;
        
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
        
        
        
        NODE = '_node';
        LINK = '_link';
        
        
        TAU = '_ta';
        NU = '_nu';
        LAMBDA = '_la';
        PI = '_pi';
        DELTA = '_de';
        
        
        QUEUE = '_queue';
        QUEUEELE = '_element';
        CONTEXT = '_context';
        
        inFluxNode = '_if';
        efFluxNode = '_ef';
        
    end

    methods 

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % basic: init,put(s),get(s),modify(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = myDS(oPath,initTable,initSequence)
            
            obj.sequence = [];
            obj.seqPtr = 0;
            
            
            if nargin >= 0
                obj.base_outPath = {''};
                obj.ds = table(zeros(0,0,'uint32'),[],{},{},zeros(0,0,'uint32'),zeros(0,0,'uint32'));
                obj.ds.Properties.VariableNames = {'key' 'time' 'data' 'type' 'linkKey' 'slinkKey'};
            end
            
            if nargin >= 1
                obj.base_outPath = oPath;
            end
            
            if nargin >= 2
                 obj.ds = initTable;
            end
            
            if nargin >= 3
                obj.sequence = initSequence;
                obj.seqPtr = 1;
            end
            
            setEmptyPtr(obj);
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put object to table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnKey,linkKey] = put(obj,data,type,linkKey,slinkKey,isSimple)
            % no source link - then set to 0
            if nargin <= 4
                slinkKey = 0;
            end
            % if 
            if nargin <= 5
                isSimple = false;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % init return key as empty
            returnKey = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % insertion logic START
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if type contains keyword image
            % then insert "self" frame for image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if contains(lower(type),'image')
                [K,linkKey] = obj.put(eye(3),'iFrame',linkKey);
                returnKey = K;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % insertion logic END
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate time-key
            [K,T] = obj.generateKeyTime();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(linkKey)
                 linkKey = K;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if slinkKey == 0
                slinkKey = K;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % perform insert
            pr = obj.emptyPtr;
            obj.ds.key(pr) = K;
            obj.ds.time(pr) = T;
            obj.ds.data{pr} = data;
            obj.ds.type{pr} = type;
            obj.ds.linkKey(pr) = linkKey;
            obj.ds.slinkKey(pr) = slinkKey;
            obj.emptyPtr = obj.emptyPtr + 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            if isempty(returnKey)
                returnKey = K;
            end
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate influx and efflux nodes
            if ~isSimple && obj.effluxRender && ~strcmp(type,[myDS.inFluxNode myDS.NODE]) && ~strcmp(type,[myDS.efFluxNode myDS.NODE])
                obj.renderEffluxNode(K);
            end
            if ~isSimple && obj.influxRender && ~strcmp(type,[myDS.inFluxNode myDS.NODE]) && ~strcmp(type,[myDS.efFluxNode myDS.NODE])
                obj.renderInfluxNode(K);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put link to table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [linkKey] = putLink(obj,sourceKey,linkType,targetKey,isDirect)
            % append key work link to link type
            linkType = [linkType myDS.LINK];
            
            
            % default is non-direct linking
            if nargin <= 4
                isDirect = [false false];
            end
            
            
            % if not direct source - then attempt to get efflux node as source
            if ~isDirect(1)
                sourceEffluxNode = obj.getEffluxNode(sourceKey);
                if ~isnan(sourceEffluxNode)
                    sourceKey = sourceEffluxNode;
                end
            end
            
            % if not direct target - then attempt to get influx node as target
            if ~isDirect(1)
                targetInfluxNode = obj.getInfluxNode(targetKey);
                if ~isnan(targetInfluxNode)
                    targetKey = targetInfluxNode;
                end
            end
            
            
            % create link
            linkKey = obj.put('',linkType,targetKey,sourceKey,true);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put node to table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnKey] = putNode(obj,data,nodeType,linkKey,isSimple)
            if nargin <= 4
                isSimple = false;
            end
            nodeType = [nodeType myDS.NODE];
            returnKey = obj.put(data,nodeType,linkKey);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [key] = generateNew(obj,type)
            switch type
                case 'nu'
                    [key] = obj.put([],myDS.NU,NaN,NaN,true);
                case 'delta'
                    [key] = obj.put([],myDS.DELTA,NaN,NaN,true);
                case 'lambda'
                    [key] = obj.put([],myDS.LAMBDA,NaN,NaN,true);
                    [ikey] = obj.put([],myDS.inFluxNode,NaN,NaN,true);
                    [ekey] = obj.put([],myDS.efFluxNode,NaN,NaN,true);
                    obj.putLink(ikey,'_lu',key,true);
                    obj.putLink(ekey,'_lu',key,true);
                case 'pi'
                    [key] = obj.put([],myDS.PI,NaN,NaN,true);
                case 'tau'
                    [tau_key] = obj.put([],myDS.TAU,NaN,NaN,true);
                    [nu_key] = obj.generateNew('nu');
                    [delta_key] = obj.generateNew('delta');
                    [lambda_key] = obj.generateNew('lambda');
                    [pi_key] = obj.generateNew('pi');
                    
                    obj.putLink(nu_key,'_fu',tau_key,true);
                    obj.putLink(delta_key,'_fu',nu_key,true);
                    obj.putLink(lambda_key,'_fu',nu_key,true);
                    obj.putLink(pi_key,'_fu',nu_key,true);
                    key = tau_key;
            end
        end
        
        function [lambdaKey] = getTauLinkComplex(obj,tauID)
            [lambdaKey] = obj.searchByTypelinkKey(tauID,'_la');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get object/link from table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [returnTable] = get(obj,keyList)
            for e = 1:numel(keyList)
                IDX(e) = find(obj.ds.key == keyList(e));
            end
            returnTable = obj.ds(IDX,:);
        end
        
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

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mass import files from list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = importImagesFromList(obj,fileList,type)
            for f = 1:numel(fileList)
                obj.put(fileList{f},type,'');
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
        % raw iget
        function [] = iget(obj,file)
        
        end
        % raw iput
        function [] = iput(obj,file)
        
        end
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
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % methods for structures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SEQUENCE START
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
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
        function [data,info] = read(obj)
            %%%%%%%%%%%%%%%%%%%%
            % get the key for the "next" object to be read
            nextKey = obj.sequence(obj.seqPtr);
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

        function [tf] = hasdata(obj)
            tf = numel(obj.sequence) > obj.seqPtr;
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % basic create and cast key <--> index
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate lambda/link node pairs for data/delta node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get flux by type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get inFlux object for node keyList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [inFluxKey] = getInfluxNode(obj,keyList)
            for e = 1:numel(keyList)
                [inFluxKey(e)] = obj.searchByTypelinkKey([myDS.inFluxNode myDS.NODE],keyList);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get inFlux object for node keyList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [efFluxKey] = getEffluxNode(obj,keyList)
            for e = 1:numel(keyList)
                [efFluxKey(e)] = obj.searchByTypelinkKey([myDS.efFluxNode  myDS.NODE],keyList);
            end
        end
        
        
        
        
        
        function [sourceKeyList] = getInfluxLinkList(obj,objKey)
            sourceKeyList = [];
            if ~obj.isType(objKey,[myDS.inFluxNode  myDS.NODE])
                sourceKeyList = obj.getFlux(objKey,1);
            end
        end
        
        function [targetKeyList] = getEffluxLinkList(obj,objKey)
            targetKeyList = [];
            if ~obj.isType(objKey,[myDS.efFluxNode  myDS.NODE])
                targetKeyList = obj.getFlux(objKey,8);
            end
        end
        
        
            
        
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