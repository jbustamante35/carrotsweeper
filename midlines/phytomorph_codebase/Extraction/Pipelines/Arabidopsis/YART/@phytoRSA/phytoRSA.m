classdef phytoRSA < handle
    properties
        pathList;
        mainRoot;
    end
    
    methods
        
        function [obj] = phytoRSA()
            
        end        
        
        function [] = addBranch(obj,branch)
            obj.pathList{end+1} = branch;
        end
        
        function [] = defineMainRoot(obj,mainRoot)
            obj.mainRoot = mainRoot;
        end
        
        function [] = moveRSA(obj,T)
            for e = 1:numel(obj.pathList)
                obj.pathList{e}.pathData = T*[obj.pathList{e}.pathData;ones(1,size(obj.pathList{e}.pathData,2))];
            end
        end
    end
end