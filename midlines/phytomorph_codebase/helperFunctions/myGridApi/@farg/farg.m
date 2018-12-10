classdef farg < variable
    properties
        argNumber;
    end
    
    properties (Constant)
         classDef = [functionObjectProps.defaultDomain functionObjectProps.defaultDomain 'farg'];
         objectDomain = [functionObjectProps.defaultDomain 'farg/'];
         
         classDefDCV = [functionObjectProps.defaultDomain 'ClassDef/directContainsValue'];
         classDefICV = [functionObjectProps.defaultDomain  'ClassDef/indirectContainsValue'];
         
         ICV_type = [functionObjectProps.defaultDomain  'indirectContainsValue/'];
         DCV_type = [functionObjectProps.defaultDomain 'directContainsValue/'];
    end
    
    methods
        function [obj] = farg(value,argNumber,valueDatabase,valueStream,toPersist)
            obj = obj@variable(value,valueDatabase,valueStream,[0 0]);
            obj.argNumber = argNumber;
            
            if toPersist
                indirectValueString = obj.persistValueToStream(value,valueStream);
                obj.persistExistanceToStream(value,valueDatabase,indirectValueString);
            end
        end
    end
    
    methods (Access = private)
        function [indirectValueString] = persistValueToStream(obj,value,valueStream)
            [indirectValueString] = persistValueToStream(value,obj.argNumber,valueStream);
        end
        function [] = persistExistanceToStream(obj,value,valueDatabase,indirectValueString)
            persistExistanceToStream(value,obj.uniqueKey,indirectValueString,valueDatabase);
        end
    end
    
    methods (Static)
        function [] = renderClass(valueDatabase)
            import com.franz.agraph.repository.AGServer
            import org.openrdf.model.vocabulary.RDF
            import org.openrdf.model.vocabulary.OWL
            import org.openrdf.model.vocabulary.RDFS
            import org.openrdf.query.QueryLanguage
            vf = valueDatabase.getRepository().getValueFactory();
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            classDef = vf.createURI(farg.classDef);
            valueDatabase.add(classDef,RDF.TYPE,RDFS.CLASS,cont);
            dcv = vf.createURI(farg.classDefDCV);
            icv = vf.createURI(farg.classDefDCV);
            valueDatabase.add(dcv,RDF.TYPE,RDFS.CLASS,cont);
            valueDatabase.add(icv,RDF.TYPE,RDFS.CLASS,cont);
            DCV_type = vf.createURI(farg.DCV_type);
            ICV_type = vf.createURI(farg.ICV_type);
            valueDatabase.add(DCV_type,RDF.TYPE,dcv,cont);
            valueDatabase.add(ICV_type,RDF.TYPE,icv,cont);
        end
    end
end