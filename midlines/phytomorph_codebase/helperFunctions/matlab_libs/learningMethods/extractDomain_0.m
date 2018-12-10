function [retDomain] = extractDomain_0(domainTypes,rangeTypes,domainGlue,domainToRangeGlue,conn)
    

    [domainTypes] = createURI(conn,domainTypes,[baseDomain.className fmo.subDomain]);
    [rangeTypes] = createURI(conn,rangeTypes,[baseDomain.className fmo.subDomain]);
    [domainGlue] = createURI(conn,domainGlue);
    [domainToRangeGlue] = createURI(conn,domainToRangeGlue);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % import java files
    import org.openrdf.model.vocabulary.RDF;
    import org.openrdf.model.vocabulary.RDFS;
    import org.openrdf.query.QueryLanguage;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % create value factory
    vf = conn.getRepository().getValueFactory();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % make null
    cont = javaArray('org.openrdf.model.Resource',1);
    cont(1) = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create a domain, range and connection-super properties
    connectionObject = vf.createURI(fo2.connectionObject);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % componets for select string
    % init select string
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    selectString = '?o1 ';
    % common select pattern with $ as numeric var 
    selectPattern = ['?o$ '];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % components for query string
    % init query string
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    queryString = '{?o1 ?typePredicate ?domainType1 . ';
    % common query pattern witb $ as numeric var and t$ has ith type var
    % o1 is glued to all objects and all objects are glued to each other
    queryPattern = ['?o1 ?domainGlue ?o$ . ?o$ ?typePredicate ?domainType$ . '];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % components for constraints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constraintPattern = '?o1 ?domainToRangeGlue ?f$ . ?f$ ?typePredicate ?rangeType$ . ';
    constraintString = '{';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build query string
    for e = 2:numel(domainTypes)
        queryString = [queryString strrep(queryPattern,'$',num2str(e))];
        selectString = [selectString strrep(selectPattern,'$',num2str(e))];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build range constraint string
    for e = 1:numel(rangeTypes)
        constraintString = [constraintString strrep(constraintPattern,'$',num2str(e))];
    end
    constraintString = [constraintString '}}'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct whole string and prepare query object
    wholeString = ['SELECT ' selectString 'WHERE ' queryString ' FILTER NOT EXISTS ' constraintString ];
    % create query object
    tupleQuery = conn.prepareTupleQuery(QueryLanguage.SPARQL, wholeString);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bind predicate(s)
    tupleQuery.setBinding('domainGlue',domainGlue);
    tupleQuery.setBinding('domainToRangeGlue',domainToRangeGlue);
    tupleQuery.setBinding('typePredicate',RDF.TYPE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bind subjects and objects
    for e = 1:numel(domainTypes)
        tupleQuery.setBinding(strrep('domainType$','$',num2str(e)),domainTypes{e});
    end
    for e = 1:numel(rangeTypes)
        tupleQuery.setBinding(strrep('rangeType$','$',num2str(e)),rangeTypes{e});
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % turn on reasoning
    tupleQuery.setIncludeInferred(logical(1));
    % evaluate query for domain object
    result = tupleQuery.evaluate();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    retDomain = {};
    while result.hasNext()
        bindingSet = result.next();
        retDomain{end+1}{1} = bindingSet.getValue('o1');
        for e = 2:numel(domainTypes)
            retDomain{end}{e} = bindingSet.getValue(strrep('o$','$',num2str(e)));
        end
    end

end

function [ret] = createURI(conn,URIstring,varargin)
    retFlag = 0;
    vf = conn.getRepository().getValueFactory();
    if nargin == 3
        URIbuilder = @(x)vf.createURI(varargin{1},x); 
    else
        URIbuilder = @(x)vf.createURI(x);
    end
    if ~iscell(URIstring)
        retFlag = 1;
        URIstring = {URIstring};
    end
    for e = 1:numel(URIstring)
        ret{e} = URIbuilder(URIstring{e});
    end
    if retFlag
        ret = ret{1};
    end
end





