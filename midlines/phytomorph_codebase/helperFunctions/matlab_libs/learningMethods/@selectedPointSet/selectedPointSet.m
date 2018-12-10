classdef selectedPointSet
    
    properties (Constant)
        % make semantics class name for reference
        className = [baseDomain.className 'matlabObjects/selectedPointSet'];
        subDomain = 'matlabObjects/';
        functionAttachment_hand = [baseDomain.className 'matlabObjects/handSelected'];
        functionAttachment_computed = [baseDomain.className 'matlabObjects/computerGenerated'];
    end
    
    properties
        pointSet = [];
        oPath = '';
        key = '';
        attachedFileURI = '';
        moniker = 'pointSet';
    end
    
    methods
        function [obj] = selectedPointSet(fileURI,varargin)
            obj.key = generateObjectKey(obj);
            obj.attachedFileURI = char(fileURI.toString());
            if nargin == 2
                obj.pointSet = varargin{1};
            end
        end
        
        function [] = setPointSet(obj,pointSet)
            obj.pointSet = pointSet;
        end
        
        function [pointSet] = getPointSet(obj)
            pointSet = obj.pointSet;
        end
        
        function [] = setFileURI(obj,fileURI)
            obj.attachedFileURI = fileURI;
        end
        
        function [fileURI] = getAttachedFileURI(obj)
            fileURI = obj.attachedFileURI;
        end
        
        function [] = persist(obj,pInfo)
            fprintf(['start persisting@' obj.moniker '@' obj.key '\n']);tic
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % create value factory
            vf = pInfo.conn.getRepository().getValueFactory();
            % persist the data to disk
            obj.oPath = pInfo.oPath;
            matFile = [pInfo.oPath obj.key '.mat'];
            save(matFile,'obj');
            % conect the object via semantic linkage
            % create image file uri
            imageFileURI = vf.createURI(obj.attachedFileURI);
            % create the pointset fileURI
            fileNameSpace = ['file:/' pInfo.oPath];
            % make a feature map object
            pointSetFile = vf.createURI(fileNameSpace,[obj.key '.mat']);
            % attach semantics
            predicate = vf.createURI(selectedPointSet.functionAttachment_hand);
            % class type semantics
            className = vf.createURI(selectedPointSet.className);
            pInfo.conn.add(pointSetFile,predicate,imageFileURI,cont);
            pInfo.conn.add(pointSetFile,RDF.TYPE,className,cont);
            
            fprintf(['end persisting@' obj.moniker '@' obj.key ':' num2str(toc) '\n']);
        end
        
        function [nKey] = generateObjectKey(obj)
            % random+clock key
            toHash = [datestr(clock) '-' num2str(randi(10000,1))];
            nKey = ['o_' num2str(string2hash(toHash))];
            fprintf(['Hashed for feature map object:' toHash '-->' nKey '\n'])
        end
    end
    
    methods (Static)
        function [] = renderClassType(conn)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            import org.openrdf.model.vocabulary.OWL;
            %% create value factory
            vf = conn.getRepository().getValueFactory();
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % make a feature map object
            fmoClass = vf.createURI(fmo.className);
            conn.add(fmoClass,RDF.TYPE,RDFS.CLASS,cont);
            % attach connection semantics to the hand slected 
            handType = vf.createURI(selectedPointSet.functionAttachment_hand);
            % create a domain, range and connection-super properties
            connectionObject = vf.createURI(fo2.connectionObject);
            % hand connection is a connection type
            conn.add(handType,RDFS.SUBPROPERTYOF,connectionObject,cont);
        end
    end


end