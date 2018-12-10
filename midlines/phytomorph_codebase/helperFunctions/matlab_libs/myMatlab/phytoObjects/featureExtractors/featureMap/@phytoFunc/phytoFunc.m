classdef phytoFunc  < matlab.mixin.Copyable
    
    properties
        para;
        notes;
        func;
    end
    
    methods
        %%%%%%%%%%%%
        % constructor
        function [obj] = phytoFunc()
            obj.para = [];
        end
        %%%%%%%%%%%%
        % get the parameters
        function [para] = getPara(obj)
            para = obj.para;
        end
        %%%%%%%%%%%%
        % set the parameters
        function [] = setPara(obj,para)
            obj.para = para;
        end
    end
    
    methods (Static)
        %%%%%%%%%%%%
        % construct crop box
        function ret = constructCropBox_atP(P,w,h)
            ret{1} = P;
            ret{2} = w;
            ret{3} = h;
        end
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute feature map at P
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fM X] = computeFeatureMap_atP(obj,X,P)
            % make call to myReader to read a crop box
            X = myReader(X,'atP',P);
            %{
            hold off
            imshow(X,[]);
            drawnow
            hold on
            %}
            % compute via self call
            fM = obj.computeFeatureMap(X);
        end
    end
        
    
    methods (Abstract)
        %%%%%%%%%%%%
        % fill in
        fM = computeFeatureMap(obj,X);       
    end
end