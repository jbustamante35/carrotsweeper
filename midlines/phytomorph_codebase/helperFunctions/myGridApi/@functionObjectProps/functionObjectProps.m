classdef functionObjectProps < handle
    properties (Constant)
        defaultDomain = 'http://www.phytomorph.wisc.edu/functionsGraph/';
        classDefLocation = 'classDef/';
    end
    
    
    properties
        uniqueKey = '';
    end
    
    methods
        function [obj] = functionObjectProps()
            obj.uniqueKey =  generateUniqueKey();
            urlID = getURLid(obj);
            valueDatabase.add(urlID, RDF.TYPE,classDef,cont);
        end
        
        function [urlID] = getURLid(obj)
            % get class
            classID = class(obj);
            % render unique key as farg object
            urlID = vf.createURI([functionObjectProps.defaultDomain classID filesep obj.uniqueKey]);
        end
    end
    
    methods (Static)
        function [classDefURL] = renderClassDef(obj,valueDatabase,OBJECT)
            if ~ischar(obj)
                obj = class(obj);
            end
            % import needed libs
            %import org.openrdf.model.vocabulary.RDF
            %import org.openrdf.model.vocabulary.RDFS
            % get the java null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % get the class der url
            classDefURL = functionObjectProps.getClassDefURLfromString(obj,valueDatabase);
            % create function object
            valueDatabase.add(classDefURL,RDF.TYPE,OBJECT,cont);
        end
        
        function [] = renderClassDefLevels(obj,valueDatabase,OBJECT,PREDICATE)
            if ~iscell(obj)
                obj = {obj};
            end
            %import org.openrdf.model.vocabulary.RDFS
            % get the java null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            for l = 1:numel(obj)
                % render the class for the object in question
                currentClassDef = functionObjectProps.renderClassDef(obj{l},valueDatabase,OBJECT);
                % list the super classes for the object
                superClassList = functionObjects.directSuperClasses(obj{l});
                % render each as a class and subclass the current object
                for e = 1:numel(superClassList)
                    % render the parent class
                    tmpClassDef = functionObjects.renderClassDef(superClassList{e},valueDatabase);
                    % render the subclass 
                    valueDatabase.add(currentClassDef,PREDICATE,tmpClassDef,cont);
                    % render the gg-parents
                    functionObjectProps.renderClassDefLevels(superClassList{e},valueDatabase,OBJECT,PREDICATE);
                end
            end
        end
        
        function [finalSuperClassList] = directSuperClasses(className)
            superClassList = superclasses(className);
            finalSuperClassList = superClassList;
            for e = 1:numel(superClassList)
                tmpList = superclasses(superClassList{e});
                finalSuperClassList = setdiff(finalSuperClassList,tmpList);
            end
        end
        
        function [classDef] = getClassDefURL(obj,valueDatabase)
            classDef = getClassDefURLfromString(class(obj),valueDatabase);
        end
        
        function [classDef] = getClassDefURLfromString(className,valueDatabase)
            classDef = [functionObjectProps.defaultDomain functionObjectProps.classDefLocation className];
            vf = valueDatabase.getRepository().getValueFactory();
            classDef = vf.createURI(classDef);
        end
        
        function [classURL] = getClassURL(obj,valueDatabase)
            classURL = [functionObjectProps.defaultDomain class(obj)];
            vf = valueDatabase.getRepository().getValueFactory();
            classURL = vf.createURI(classURL);
        end
    end
end