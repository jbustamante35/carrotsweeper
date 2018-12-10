classdef phytoGeoSet < myHS_X & geometrical
    % construct for geo-objects
    properties
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoGeoSet(varargin)
            % constructor
            obj@myHS_X();
            obj@geometrical();
            % if one then should be element type
            if nargin >= 1
                % set the element type
                obj.addAllow(varargin{1});
                % type operator for extraction
                %obj.extractFunction = typeOp(obj.eType,class(obj));
            end
        end
    end
end