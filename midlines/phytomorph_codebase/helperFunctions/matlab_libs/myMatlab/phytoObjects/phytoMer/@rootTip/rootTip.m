classdef rootTip < matlab.mixin.Copyable
    
    properties
        % anchor points
        tip_point;
        base_point;
        upper_base_point;
        lower_base_point
        % curves
        contour;
        midline_1;
        midline_2;
        % tip reference frame
        tip_reference_frame;
        % associated image
        image;
        flipD;
    end
    
    methods
        % constructor
        function [obj] = rootTip(varargin)
            % points
            obj.base_point = phytoApoint();
            obj.tip_point = phytoApoint();
            obj.upper_base_point = phytoApoint();
            obj.lower_base_point = phytoApoint();
            % curves
            obj.contour = phytoAcurve();
            obj.midline_1 = phytoAcurve();
            obj.midline_2 = phytoAcurve();
            % tip reference frame
            obj.tip_reference_frame = phytoAaffine();
            % image
            obj.image = imageFile();
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set function
        function [] = set(obj,fieldName,data)
            obj.(fieldName) = data;
        end
        % set data function
        function [] = setData(obj,fieldName,data)
            obj.(fieldName).setData(data);
        end
        % get function
        function [data] = get(obj,fieldName)
            data = obj.(fieldName);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set image
        function [] = setImage(obj,imageFile)
            obj.image = imageFile;
        end
        % set image data
        function [] = setImageName(obj,imageFileData)
            obj.image.setFileName(imageFileData);
        end
        % get image
        function [imageFile] = getImage(obj)
            imageFile = obj.image;
        end
        
    end
    
    

end