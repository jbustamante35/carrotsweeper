classdef funcObject 
    properties
        
    end
    
    properties (Constant)
         classDef = [functionObjects.defaultDomain 'ClassDef/funcObject'];
         classDefContain = [functionObjects.defaultDomain 'ClassDef/isInputTo'];
    end
    
    methods
        function [obj] = fargList(values,valueDatabase,valueStream)
            for e = 1:numel(values)
                farg(values{e},e,valueDatabase,valueStream,1);
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
            
            % create function object
            classDef = vf.createURI(fargList.classDef);
            valueDatabase.add(classDef,RDF.TYPE,RDFS.CLASS,cont);
            
            % create isInputTo connection type
            classDefContains = vf.createURI(fargList.classDefContains);
            valueDatabase.add(classDefContains,RDF.TYPE,RDFS.CLASS,cont);
            
        end
    end
end