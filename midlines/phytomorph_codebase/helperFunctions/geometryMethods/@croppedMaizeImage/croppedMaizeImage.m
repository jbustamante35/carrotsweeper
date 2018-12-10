classdef croppedMaizeImage < maizeImage
    properties
        cropBox;
    end
    
    methods
        function [obj] = croppedMaizeImage(filename)
            obj@maizeImage(filename);
        end
    end
end