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
%% 
InArgumentList = vf.createURI(psNameSpace,'inputArgumentList');
conn.add(InArgumentList,RDFS.SUBCLASSOF,logicalContainer,cont);

OutArgumentList = vf.createURI(psNameSpace,'ouputArgumentList');
conn.add(OutArgumentList,RDFS.SUBCLASSOF,logicalContainer,cont);


input_argument = vf.createURI(psNameSpace,'inputArgument');
conn.add(input_argument,RDFS.SUBCLASSOF,object,cont);

output_argument = vf.createURI(psNameSpace,'ouputArgument');
conn.add(output_argument,RDFS.SUBCLASSOF,object,cont);


functionObject = vf.createURI(psNameSpace,'function');

isInputTo = vf.createURI(psNameSpace,'isInputTo');
conn.add(isInputTo,RDF.TYPE,RDF.PREDICATE,cont);
conn.add(isInputTo,RDF.TYPE,RDFS.CLASS,cont);


specialOne = vf.createURI(psNameSpace,'specialOne');
conn.add(specialOne,RDFS.SUBCLASSOF,isInputTo,cont);
conn.add(specialOne,RDFS.SUBCLASSOF,container,cont);


%%



ra1 = vf.createURI(psNameSpace,'ra1');
conn.add(ra1,RDF.TYPE,isInputTo,cont);


a1 = vf.createURI(psNameSpace,'a1');
conn.add(a1,RDF.TYPE,input_argument,cont);
a2 = vf.createURI(psNameSpace,'a2');
conn.add(a2,RDF.TYPE,input_argument,cont);


arg1 = vf.createURI(psNameSpace,'arg1');
conn.add(a1,RDF.TYPE,specialOne,cont);
arg2 = vf.createURI(psNameSpace,'arg2');
conn.add(a2,RDF.TYPE,specialOne,cont);


ra1 = vf.createURI(psNameSpace,'ra1');
conn.add(ra1,RDF.TYPE,isInputTo,cont);

func = vf.createURI(psNameSpace,'func');
inSet = vf.createURI(psNameSpace,'inSet');
conn.add(inSet,isInputTo,func,cont);
%%
%% simple test

tt = vf.createURI(psNameSpace,'tt');
conn.add(tt,RDFS.SUBCLASSOF,RDF.TYPE,cont);
conn.add(tt,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);

pred = vf.createURI(psNameSpace,'pred');
conn.add(pred,RDF.TYPE,RDF.PREDICATE,cont);
conn.add(pred,RDF.TYPE,OWL.TRANSITIVEPROPERTY,cont);
conn.add(isInputTo,RDF.TYPE,RDFS.CLASS,cont);

rpred1 = vf.createURI(psNameSpace,'rpred1');
conn.add(rpred1,tt,pred,cont);
rpred2 = vf.createURI(psNameSpace,'rpred2');
conn.add(rpred2,tt,pred,cont);


o1 = vf.createURI(psNameSpace,'o1');
o2 = vf.createURI(psNameSpace,'o2');
o3 = vf.createURI(psNameSpace,'o3');
conn.add(o1,rpred1,o2,cont);
conn.add(o2,rpred2,o3,cont);



%% search facts
import org.openrdf.query.QueryLanguage
queryString = 'SELECT ?s ?p ?o   WHERE {?s ?p ?o .}';
tupleQuery = conn.prepareTupleQuery(QueryLanguage.SPARQL, queryString);
tupleQuery.setIncludeInferred(logical(1));
tupleQuery.setBinding('s',o1);
result = tupleQuery.evaluate();
while result.hasNext()
    bindingSet = result.next();
    s = bindingSet.getValue('s');
    p = bindingSet.getValue('p');
    o = bindingSet.getValue('o');
    [char(s.toString()) '-->' char(p.toString()) '-->' char(o.toString())]
end


%% make abstract production system
% create year
year = vf.createURI(psNameSpace,'year');
conn.add(year,RDFS.SUBCLASSOF,logicalContainer,cont);
% create season
season = vf.createURI(psNameSpace,'season');
conn.add(season,RDFS.SUBCLASSOF,logicalContainer,cont);
y_contains_s = vf.createURI(psNameSpace,['contains_year_season']);
conn.add(year,y_contains_s,season,cont);
conn.add(year,OWL.OBJECTPROPERTY,y_contains_s);
conn.add(y_contains_s,RDFS.SUBPROPERTYOF,contains,cont);
% create site
site = vf.createURI(psNameSpace,'site');
conn.add(site,RDFS.SUBCLASSOF,logicalContainer,cont);
s_contains_si = vf.createURI(psNameSpace,['contains_season_site']);
conn.add(season,s_contains_si,site,cont);
conn.add(season,OWL.OBJECTPROPERTY,s_contains_si);
conn.add(s_contains_si,RDFS.SUBPROPERTYOF,contains,cont);
% create field
field = vf.createURI(psNameSpace,'field');
conn.add(field,RDFS.SUBCLASSOF,logicalContainer,cont);
si_contains_f = vf.createURI(psNameSpace,['contains_site_field']);
conn.add(site,si_contains_f,field,cont);
conn.add(site,OWL.OBJECTPROPERTY,si_contains_f);
conn.add(si_contains_f,RDFS.SUBPROPERTYOF,contains,cont);
% create range
range = vf.createURI(psNameSpace,'range');
conn.add(range,RDFS.SUBCLASSOF,logicalContainer,cont);
f_contains_r = vf.createURI(psNameSpace,['contains_field_range']);
conn.add(field,f_contains_r,range,cont);
conn.add(field,OWL.OBJECTPROPERTY,f_contains_r);
conn.add(f_contains_r,RDFS.SUBPROPERTYOF,contains,cont);
% create plot
plot = vf.createURI(psNameSpace,'plot');
conn.add(plot,RDFS.SUBCLASSOF,logicalContainer,cont);
r_contains_p = vf.createURI(psNameSpace,['contains_range_plot']);
conn.add(range,r_contains_p,plot,cont);
conn.add(range,OWL.OBJECTPROPERTY,r_contains_p);
conn.add(r_contains_p,RDFS.SUBPROPERTYOF,contains,cont);
% create row
row = vf.createURI(psNameSpace,'row');
conn.add(row,RDFS.SUBCLASSOF,logicalContainer,cont);
p_contains_r = vf.createURI(psNameSpace,['contains_plot_row']);
conn.add(plot,p_contains_r,row,cont);
conn.add(plot,OWL.OBJECTPROPERTY,p_contains_r);
conn.add(p_contains_r,RDFS.SUBPROPERTYOF,contains,cont);