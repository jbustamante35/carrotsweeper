%% import the libraries
import com.franz.agraph.repository.AGServer
import org.openrdf.model.vocabulary.RDF
import org.openrdf.model.vocabulary.OWL
import org.openrdf.model.vocabulary.RDFS
import org.openrdf.query.QueryLanguage
%% connection data
SERVER_URL = 'http://localhost:10035';
CATALOG_ID = 'pipelines';
REPOSITORY_ID = 'functionCalls';
USERNAME = 'devUser';
PASSWORD = 'devUser';
%% connect
server = AGServer(SERVER_URL, USERNAME, PASSWORD);
server.listCatalogs()
catalog = server.getRootCatalog();
%% create catalog
myRepository = catalog.createRepository(REPOSITORY_ID);
%% init
myRepository.initialize();
myRepository.isWritable();
%% get connection
conn = myRepository.getConnection();
myRepository.getRepositoryID()
%% create value factory
vf = conn.getRepository().getValueFactory();
%% make null
cont = javaArray('org.openrdf.model.Resource',1);
cont(1) = [];
%% make grand name space
psNameSpace = 'http://example.org/pipeLines/';
%% make container system
% make super property the inverse of sub property
super_prop = vf.createURI(psNameSpace,'SUPER_PROPERTYOF');
conn.add(super_prop, OWL.INVERSEOF, RDFS.SUBPROPERTYOF,cont);

% make a container which is-a object
container = vf.createURI(psNameSpace,'container');
object = vf.createURI(psNameSpace,'object');
conn.add(container,RDFS.SUBCLASSOF,object,cont);
conn.add(container,RDF.TYPE,RDF.SUBJECT,cont);
conn.add(container,RDF.TYPE,RDFS.CLASS,cont);


% make a transitive inference
tInfers = vf.createURI(psNameSpace,'transitive_infers');
conn.add(tInfers,RDF.TYPE,RDF.PREDICATE,cont);
conn.add(tInfers,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);


% make contains
contains = vf.createURI(psNameSpace,'contains');
i_contain = vf.createURI(psNameSpace,'containedIn');
conn.add(contains,RDF.TYPE,RDF.PREDICATE,cont);
conn.add(i_contain, OWL.INVERSEOF, contains,cont);
conn.add(contains,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
%conn.add(contains,RDFS.DOMAIN,container,cont);
%conn.add(contains,RDFS.RANGE,object,cont);



% make a hasa
hasa = vf.createURI(psNameSpace,'has-a');
conn.add(hasa,RDF.TYPE,RDF.PREDICATE,cont);
conn.add(contains,RDF.TYPE,RDF.PROPERTY,cont);


% make logical and physical containers
logicalContainer = vf.createURI(psNameSpace,'logicalContainer');
physicalContainer = vf.createURI(psNameSpace,'physicalContainer');
conn.add(logicalContainer,RDFS.SUBCLASSOF,container,cont);
conn.add(physicalContainer,RDFS.SUBCLASSOF,container,cont);
%% simple static property distribution
p1 = vf.createURI(psNameSpace,'p1');
conn.add(p1,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
p2 = vf.createURI(psNameSpace,'p2');
conn.add(p2,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
p3 = vf.createURI(psNameSpace,'p3');
conn.add(p3,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
% create values
v1 = vf.createURI(psNameSpace,'v1');
v2 = vf.createURI(psNameSpace,'v2');
v3 = vf.createURI(psNameSpace,'v3');
% create master property staticly
mp = vf.createURI(psNameSpace,'mp');
conn.add(mp,RDFS.SUBPROPERTYOF,p1,cont);
conn.add(mp,RDFS.SUBPROPERTYOF,p2,cont);
conn.add(mp,RDFS.SUBPROPERTYOF,p3,cont);
% create property list
pL = vf.createURI(psNameSpace,'propertyList');
conn.add(pL,p1,v1,cont);
conn.add(pL,p2,v2,cont);
conn.add(pL,p3,v3,cont);
testObject = vf.createURI(psNameSpace,'testObject');
% apply property list to test object and then test for distribution
conn.add(testObject,mp,pL,cont);
%% dynamic property distribtuion
p1 = vf.createURI(psNameSpace,'p1');
conn.add(p1,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
p2 = vf.createURI(psNameSpace,'p2');
conn.add(p2,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
p3 = vf.createURI(psNameSpace,'p3');
conn.add(p3,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);


BP = vf.createURI(psNameSpace,'bulkProperty');
mp = vf.createURI(psNameSpace,'mp');


conn.add(BP,RDFS.SUBPROPERTYOF,p2,cont);
conn.add(BP,RDFS.SUBPROPERTYOF,p3,cont);
conn.add(mp,RDFS.SUBPROPERTYOF,BP,cont);


% create property list
pL = vf.createURI(psNameSpace,'propertyList');
conn.add(pL,p1,v1,cont);
conn.add(pL,p2,v2,cont);
conn.add(pL,p3,v3,cont);
testObject = vf.createURI(psNameSpace,'testObject');
% apply property list to test object and then test for distribution
conn.add(testObject,mp,pL,cont);


fprintf(['TEST START property2\n']);
result = conn.getStatements(testObject, [], v2 , logical(1),cont);
while result.hasNext()
    result.next()
end

fprintf(['TEST START property1\n']);
result = conn.getStatements(testObject, [], v1 , logical(1),cont);
while result.hasNext()
    result.next()
end

conn.add(BP,RDFS.SUBPROPERTYOF,p1,cont);

fprintf(['TEST START property2\n']);
result = conn.getStatements(testObject, [], v2 , logical(1),cont);
while result.hasNext()
    result.next()
end

fprintf(['TEST START property1\n']);
result = conn.getStatements(testObject, [], v1 , logical(1),cont);
while result.hasNext()
    result.next()
end

%% 
fprintf(['TEST START\n']);
import org.openrdf.query.QueryLanguage
queryString = 'SELECT ?s ?p ?o   WHERE {?s ?p ?o .}';
tupleQuery = conn.prepareTupleQuery(QueryLanguage.SPARQL, queryString);
tupleQuery.setIncludeInferred(logical(1));
%tupleQuery.setBinding('s',p1);
%tupleQuery.setBinding('p',bulkProp);
%tupleQuery.setBinding('p',RDF.TYPE);
%tupleQuery.setBinding('p',contains);
result = tupleQuery.evaluate();
while result.hasNext()
    bindingSet = result.next();
    s = bindingSet.getValue('s');
    p = bindingSet.getValue('p');
    o = bindingSet.getValue('o');
    [char(s.toString()) '-->' char(p.toString()) '-->' char(o.toString())]
end
%%

fprintf(['TEST START\n']);
result = conn.getStatements([], bulkProp, [] , logical(1),cont);
while result.hasNext()
    result.next()
end







