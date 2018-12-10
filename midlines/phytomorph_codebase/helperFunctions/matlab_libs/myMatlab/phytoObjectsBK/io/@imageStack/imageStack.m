classdef imageStack < myHS_X
    properties
        
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = imageStack(varargin)
            obj@myHS_X('imageFile');
            if nargin == 1
                for e = 1:numel(varargin{1})
                    iF = imageFile(varargin{1}{e});
                    putImage(obj,iF);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get image from list
        function [iF] = getImage(obj,n)
            iF = getElement(obj,n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put image into list
        function [] = putImage(obj,iF,n)
            if nargin == 2
                n = numel(obj.S)+1;
            end
            putElement(obj,iF,n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sort
        function [] = sort(obj)
            try
                for e = 1:numel(obj)
                    cur = getElement(obj,e);
                    nm(e) = str2num(cur.getName()); 
                end
                [j sidx] = sort(nm);
            catch ME
                
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % import from cell array
    methods (Static)
        function [isS] = imageStackSet(set)
            isS = myHS_X('imageStack');
            for e = 1:numel(set)
                isS{e} = imageStack(set{e});
            end
            
        end
    end
end