classdef deformationMap  < phytoFunc
    
    properties
        
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = deformationMap(varargin)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init func
                %%%%%%%%%%%%%%%%%%%%%%%%%%%    
            obj.func = @(image,para)miniTrack(image,para);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init notes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.notes = 'function call for generating deformation about P';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init default values for para
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % I(:,:,1)      is the first image
            % I(:,:,2)      is the second image
            % para.P        is the point to track
            % para.domain   is the tracking domain
            % para.RADIUS   is the clipping window
            % para.init_T   is the init vector for transformation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                
            obj.para.domain.value = [];
            obj.para.vars.init_T = [];           
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init non-default parameters
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin >= 1;  obj.para.vars.sig.value = varargin{1};end;            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute feature map
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fM] = computeFeatureMap(obj,X)
            fM = obj.func(X,obj.para);
        end
    end
end