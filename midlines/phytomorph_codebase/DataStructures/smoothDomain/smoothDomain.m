classdef smoothDomain < smoothTensor

    properties
    end
    
    methods
        function [obj] = smoothDomain(varargin)
            obj = obj@smoothTensor(varargin{:});
        end
        
    end


end