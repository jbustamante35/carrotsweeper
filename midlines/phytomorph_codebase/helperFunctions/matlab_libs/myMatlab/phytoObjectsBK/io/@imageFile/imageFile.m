classdef imageFile < matlab.mixin.Copyable
    
    properties
        fileName;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = imageFile(varargin)
            obj.fileName = '';
            if nargin == 1
                obj.fileName = varargin{1};
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set file name
        function [] = setFileName(obj,fileName)
            obj.fileName = fileName;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get file name
        function [fileName] = getFileName(obj)
            fileName = obj.fileName;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read
        function [I] = read(obj)
            I = myReader(obj.fileName);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get file name
        function [nm] = getName(obj)
            [pth,nm,ext] = fileparts(obj.fileName);
        end
    end
end