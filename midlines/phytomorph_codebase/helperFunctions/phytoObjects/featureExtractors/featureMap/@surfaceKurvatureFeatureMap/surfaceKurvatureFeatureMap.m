classdef surfaceKurvatureFeatureMap  < phytoFunc
    
    properties
        
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = surfaceKurvatureFeatureMap(varargin)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init func
                %%%%%%%%%%%%%%%%%%%%%%%%%%%    
            obj.func = @(image,para)surKur(image,para); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init notes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.notes = 'function call for generating curvature for image as surface';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init values for para
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.para.scales.value = 1;
            obj.para.scales.notes = 'scales for computing curvature';
            obj.para.resize.value = 1;
            obj.para.resize.notes = 'resize amount for image';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init values for para
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin >= 1;  obj.para.scales.value = varargin{1};end;
            if nargin >= 2;  obj.para.resize.value = varargin{2};end;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute feature map
        %%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function [fM] = computeFeatureMap(obj,X)
            fM = obj.func(X,obj.para);
        end        
    end
end