classdef fo2 < matlab.mixin.Copyable%handle
    % function object
    properties
        % 
        key;
        moniker;            % human name of function
        domainType;         % NOT NEEDED - replaced by linking function
        rangeType;          % NOT NEEDED - replaced by linking function
        toRepresent;        % NOT NEEDED?
        
        domainExtract_func; % x = domainExtract_func(conn);
        func_handle;        % y = f(x);
        linking_func;       % l(y) <-- y = f(x)
        
        
        semanticType;       
    end
    
     properties (Constant)
        % make semantics class name for reference
        className = [baseDomain.className 'matlabObjects/functionObject'];
        % make semantics class name for reference
        exeType = [baseDomain.className 'matlabObjects/functionObject_exeType'];
        % make base for exe key
        exeBase = [baseDomain.className 'matlabObjects/functionObject_exe_'];
        domain = [baseDomain.className 'matlabObjects/domain'];
        domainBase = [baseDomain.className 'matlabObjects/domain_'];
        domainSize = [baseDomain.className 'matlabObjects/domain_size'];
        rangeBase = [baseDomain.className 'matlabObjects/range_'];
        range = [baseDomain.className 'matlabObjects/range'];
        rangeSize = [baseDomain.className 'matlabObjects/range_size'];
        connectionObject = [baseDomain.className 'matlabObjects/connected'];
    end
    
    methods
        function [obj] = fo2(domainExtract_func,moniker,func_handle,linking_func)
            obj.domainExtract_func = domainExtract_func;
            %obj.domainType = domainType;
            obj.moniker = moniker;
            obj.func_handle = func_handle;
            %obj.rangeType = rangeType;
            obj.linking_func = linking_func;
            %obj.toRepresent = toRepresent;
        end
        
        function [x] = reifyDomain(obj,conn)
            x = obj.domainFunc(conn);
        end
        
        function [varargout] = run(obj,pInfo,varargin)
           
            
          
            
            % run the function
            varargout = obj.func_handle(varargin{:});
            
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            
            
            
            % create value factory
            vf = pInfo.conn.getRepository().getValueFactory();
            % make a feature map object
            foType = vf.createURI(obj.semanticType);
            % make a executable object
            %exeType = vf.createURI(fo2.exeType);
            % make a executable object
            exeInstance = vf.createURI(fo2.exeBase,fo2.getExeKey());
            
            
            % register the connection from domain{n}->func->codomain{m}
            pInfo.conn.add(exeInstance,exeType,foType,cont);
            
            % loop over outputs and atach to exeInstance
            for e = 1:numel(varargout)
                % addatch moniker to the featureMapObject
                varargout{e}.setMoniker(obj.rangeType{e});
                % render the file object to disk
                codomainObject = varargout{e}.persist(pInfo);
                % register the object as fmo subtype
                fmo.renderSubClassObject(pInfo.conn,obj.rangeType{e},codomainObject);
                % make the eth range predicate
                rangePredicate = vf.createURI([fo2.rangeBase num2str(e)]);
                % register the connection from domain{n}->func->codomain{m}
                pInfo.conn.add(exeInstance,rangePredicate,codomainObject,cont);
            end
            
            
            
            % loop over inputs and attach to exeInstance
            for e = 1:numel(varargin)
                 % make the eth domain predicate
                domainPredicate = vf.createURI([fo2.domainBase num2str(e)]);
                % register the connection from domain{n}->func->codomain{m}
                pInfo.conn.add(varargin{e},domainPredicate,exeInstance,cont);
            end
        end
        
        function [varargout] = subsref(obj,s)
            flag = 1;
            if flag
                if nargout == 0
                    builtin('subsref',obj,s);
                else
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                     %= {varargout};
                end
            end
            
            %{
            flag = 1;
            
            if all(s(1).type == '()')
                if numel(s(1).subs) == 1
                    if isa(s(1).subs{1},'fmo')
                        flag = 0;
                        for e = 1:nargout
                            varargout{e} = fmo(obj.ptDetector,obj.fmo_output_moniker{e},'type');
                            obj.fmo_output_type{e} = varargout{e}.typeKey;
                        end
                    end
                end         
            end
            
            if flag
                if nargout == 0
                    builtin('subsref',obj,s);
                else
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                     %= {varargout};
                end
            end
            %}
        end
        
        function [] = persist(obj,conn,oPath)
            % import the jars
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % create value factory
            vf = conn.getRepository().getValueFactory();
            fileNameSpace = ['file:/' oPath];
            % make a feature map object
            foType = vf.createURI(fileNameSpace,[obj.moniker '.mat']);
            % register the function
            conn.add(foType,RDF.TYPE,vf.createURI(fo2.className),cont);
            
            %{
            % make the domain type and register the domain size
            domainConnection = vf.createURI(fo2.domain);
            % make domain size
            domainSize = vf.createURI(fo2.domainSize);
            dSZ = vf.createLiteral(num2str(numel(obj.domainType)));
            conn.add(foType,domainSize,dSZ,cont);
            for e = 1:numel(obj.domainType)
                domain = vf.createURI(baseDomain.className,['matlabObjects/' obj.domainType{e}]);
                % make the eth domain predicate
                domainPredicate = vf.createURI([fo2.domainBase num2str(e)]);
                conn.add(domainPredicate,RDFS.SUBPROPERTYOF,domainConnection,cont);
                conn.add(foType,domainPredicate,domain,cont);
            end
            %}
            
            %{
            % make the range type
            rangeConnection = vf.createURI(fo2.range);
            rangeSize = vf.createURI(fo2.rangeSize);
            rSZ = vf.createLiteral(num2str(numel(obj.rangeType)));
            conn.add(foType,rangeSize,rSZ,cont);
            for e = 1:numel(obj.rangeType)
                % make range object(s)
                range = vf.createURI(baseDomain.className,['matlabObjects/' obj.rangeType{e}]);
                % make range predicate
                rangePredicate = vf.createURI([fo2.rangeBase num2str(e)]);
                conn.add(rangePredicate,RDFS.SUBPROPERTYOF,rangeConnection,cont);
                conn.add(foType,rangePredicate,range,cont);
                % register the object type
                fmo.renderSubClassType(conn,obj.rangeType{e},obj.toRepresent(e));
            end
            %}
            
            
            
            obj.semanticType = char(foType.toString());
            %connectionObject = vf.createURI(fo2.connectionObject);
            %conn.add(foType,RDFS.SUBPROPERTYOF,connectionObject,cont);
            save([oPath obj.moniker '.mat'],'obj');
        end
        
        function [n] = getNargout(obj)
            n = numel(obj.fmo_output_type{e});
        end
        
        function [t] = getSemanticType(obj)
        
        end
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % renderClassType: render the needed semantics for the class type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = renderClassType(conn)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            import org.openrdf.model.vocabulary.OWL;
            % create value factory
            vf = conn.getRepository().getValueFactory();
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % create a function class object
            foClass = vf.createURI(fo2.className);
            % add class objects
            conn.add(foClass,RDF.TYPE,RDFS.CLASS,cont);
            % create a domain, range and connection-super properties
            connectionObject = vf.createURI(fo2.connectionObject);
            domainConnection = vf.createURI(fo2.domain);
            rangeConnection = vf.createURI(fo2.range);
            conn.add(connectionObject,OWL.INVERSEOF,connectionObject,cont);
            conn.add(connectionObject,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
            conn.add(domainConnection,RDFS.SUBPROPERTYOF,connectionObject,cont);
            conn.add(rangeConnection,RDFS.SUBPROPERTYOF,connectionObject,cont);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loadAndRun: run the function (f) with the domain (x) --> f(x)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = loadRunLink(f,x,pInfo)
            f = fo2.load(f);
            f.run(pInfo,x{:});
        end
        
        function [f] = load(f)
            if ~ischar(f)
                f = char(f.toString());
            end
            f = f(6:end);
            f = load(f);
            f = f.obj;
        end
        
        
        function [domainType] = getDomainType(conn,foObject,domainSize)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            import org.openrdf.query.QueryLanguage
            % select and remove those which are operated on
            vf = conn.getRepository().getValueFactory();
            queryString = 'SELECT ?o WHERE {?s ?p ?o.}';
            tupleQuery = conn.prepareTupleQuery(QueryLanguage.SPARQL, queryString);
            tupleQuery.setBinding('s',foObject);
            domainSize = char(domainSize.toString());
            domainSize(1) = [];
            domainSize(end) = [];
            domainSize = str2num(domainSize);
            for e = 1:domainSize
                domainPredicate = vf.createURI([fo2.domainBase num2str(e)]);
                tupleQuery.setBinding('p',domainPredicate);
                result = tupleQuery.evaluate();
                bindingSet = result.next();
                domainType{e}  = bindingSet.getValue('o');
            end
        end
        
        function [rangeType] = getRangeType(conn,foObject,rangeSize)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            import org.openrdf.query.QueryLanguage
            % select and remove those which are operated on
            vf = conn.getRepository().getValueFactory();
            queryString = 'SELECT ?o WHERE {?s ?p ?o.}';
            tupleQuery = conn.prepareTupleQuery(QueryLanguage.SPARQL, queryString);
            tupleQuery.setBinding('s',foObject);
            rangeSize = char(rangeSize.toString());
            rangeSize(1) = [];
            rangeSize(end) = [];
            rangeSize = str2num(rangeSize);
            for e = 1:rangeSize
                rangePredicate = vf.createURI([fo2.rangeBase num2str(e)]);
                tupleQuery.setBinding('p',rangePredicate);
                result = tupleQuery.evaluate();
                bindingSet = result.next();
                rangeType{e}  = bindingSet.getValue('o');
            end
        end
        
        function [retDomain] = getConnectedDomains(conn,domainTypes,rangeTypes)
            % import java files
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            import org.openrdf.query.QueryLanguage;
            % create value factory
            vf = conn.getRepository().getValueFactory();
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % create a domain, range and connection-super properties
            connectionObject = vf.createURI(fo2.connectionObject);
            queryString = '{?o1 ?p1 ?t1 . ';
            selectString = '?o1 ';
            queryPattern = ['?o1 ?p ?o$ . ?o$ ?p1 ?t$ . '];
            selectPattern = ['?o$ '];
            % range string
            constraintPattern = '?o1 ?p ?f$ . ?f$ ?p1 ?r$ . ';
            constraintString = '{';
            % build query string
            for e = 2:numel(domainTypes)
                queryString = [queryString strrep(queryPattern,'$',num2str(e))];
                selectString = [selectString strrep(selectPattern,'$',num2str(e))];
            end
            % build range constraint string
            for e = 1:numel(rangeTypes)
                constraintString = [constraintString strrep(constraintPattern,'$',num2str(e))];
            end
            constraintString = [constraintString '}}'];
            %queryString = [queryString '}'];
            wholeString = ['SELECT ' selectString 'WHERE ' queryString ' FILTER NOT EXISTS ' constraintString ];
            % create query object
            tupleQuery = conn.prepareTupleQuery(QueryLanguage.SPARQL, wholeString);
            % bind objects across query object
            tupleQuery.setBinding('p',connectionObject);
            tupleQuery.setBinding('p1',RDF.TYPE);
            for e = 1:numel(domainTypes)
                tupleQuery.setBinding(strrep('t$','$',num2str(e)),domainTypes{e});
            end
            for e = 1:numel(rangeTypes)
                tupleQuery.setBinding(strrep('r$','$',num2str(e)),rangeTypes{e});
            end
            % turn on reasoning
            tupleQuery.setIncludeInferred(logical(1));
            % evaluate query for domain object
            result = tupleQuery.evaluate();
            
            retDomain = {};
            while result.hasNext()
                bindingSet = result.next();
                retDomain{end+1}{1} = bindingSet.getValue('o1');
                for e = 2:numel(domainTypes)
                    retDomain{end}{e} = bindingSet.getValue(strrep('o$','$',num2str(e)));
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getExeKey: get execution key
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [exeKey] = getExeKey()
            % random+clock key
            toHash = [datestr(clock) '-' num2str(randi(10000,1))];
            exeKey = [num2str(string2hash(toHash))];
            fprintf(['Hashed for exeInstance object:' toHash '-->' exeKey '\n'])
        end
    end
end