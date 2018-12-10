classdef fargList < functionObjects
    properties
        
    end
    
    properties (Constant)
         classDef = [functionObjectProps.defaultDomain 'ClassDef/fargList'];
         classDefContain = [functionObjectProps.defaultDomain 'ClassDef/containsArgument'];
         contain_type = [functionObjectProps.defaultDomain 'type/containsArgument'];
    end
    
    methods
        function [obj] = fargList(values,valueDatabase,valueStream)
            obj = obj@functionObjects();
            
            contains = vf.createURI(farg.classDefContain);
            for e = 1:numel(values)
                tmpObj = farg(values{e},e,valueDatabase,valueStream,1);
                objectID = vf.createURI(farg.objectDomain,tmpObj.uniqueKey);
                
            end 
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
            % create class def forobject
            classDef = vf.createURI(fargList.classDef);
            valueDatabase.add(classDef,RDF.TYPE,RDFS.CLASS,cont);
            % create contains class
            classDefContains = vf.createURI(fargList.classDefContain);
            valueDatabase.add(classDefContains,RDF.TYPE,RDFS.CLASS,cont);
            % create contains type
            contains = vf.createURI(fargList.contain_type);
            valueDatabase.add(contains,RDF.TYPE,classDefContains,cont);
        end
    end
end