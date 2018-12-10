classdef pImage < handle

    
    properties
        image = [];
        cache = false;
        fileLocation;
    end
    
    methods
        function [obj] = pImage(fileLocation)
            obj.fileLocation = fileLocation;
        end
        
        function [d] = readData(obj,index)
            if isempty(obj.image)
                obj.image = imread(obj.fileLocation);
            end
            
            d = obj.image(index);
        end
        
        
    end
end

%{
    for e = 1:1000
        fL{e} = FileList{1};
    end
    imds = imageDatastore(fL);


%}