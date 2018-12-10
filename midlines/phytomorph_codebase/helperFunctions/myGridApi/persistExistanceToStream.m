function [] = persistExistanceToStream(value,uniqueKey,indirectValueString,valueDatabase)
    import com.franz.agraph.repository.AGServer
    import org.openrdf.model.vocabulary.RDF
    import org.openrdf.model.vocabulary.OWL
    import org.openrdf.model.vocabulary.RDFS
    import org.openrdf.query.QueryLanguage
    vf = valueDatabase.getRepository().getValueFactory();
    cont = javaArray('org.openrdf.model.Resource',1);
    cont(1) = [];
    
    
    % render unique key as farg object
    classDef = vf.createURI(farg.classDef);
    object = vf.createURI(farg.objectDomain,uniqueKey);
    valueDatabase.add(object, RDF.TYPE,classDef,cont);
    
    % render indirect contains
    indirectContainment = vf.createURI(farg.ICV_type);
    indirectValueString = vf.createURI(indirectValueString);
    valueDatabase.add(object, indirectContainment,indirectValueString ,cont);
    
    % if char then render directly to database
    if ischar(value)
        directContainment = vf.createURI(farg.DCV_type);
        valueURL = vf.createURI(value);
        valueDatabase.add(object,directContainment,valueURL);
    end
    
end
