classdef pointExtractor  < matlab.mixin.Copyable
    
    properties
        para;
        funcList;
        notes;
    end
    
    methods
        
        % constructor
        function [obj] = pointExtractor()
            obj.para = [];
        end
        % get the parameters
        function [para] = getPara(obj)
            para = obj.para;
        end
        % set the parameters
        function [] = setPara(obj,para)
            obj.para = para;
        end
            
    end
    
    methods (Abstract)
        pL = extractPoints(obj,X);
    end
end