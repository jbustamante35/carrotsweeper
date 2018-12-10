classdef geometry2D < smoothDomain

    properties
    end
    
    methods
        function [obj] = geometry2D(varargin)
            obj = obj@smoothDomain(varargin{:});
        end
        
        
        function [] = plot(obj,CL)
            if nargin == 1
                CL = {'r' 'r'};
            end
            skip = [3 3];
            for e1 = 1:skip(1):size(obj.M,2)
                plot(obj.M(:,e1,1),obj.M(:,e1,2),CL{1})
                hold on
            end
            plot(obj.M(:,end,1),obj.M(:,end,2),CL{1})
            
            
            skip = [3 3];
            for e1 = 1:skip(1):size(obj.M,1)
                plot(obj.M(e1,:,1),obj.M(e1,:,2),CL{1})
                hold on
            end
            plot(obj.M(end,:,1),obj.M(end,:,2),CL{1})
            
            
            
            
        end
    end


end