classdef splineFit < modelFit
    
    properties
        knots;
        degree;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = splineFit(knots,degree)
            obj.knots = knots;
            obj.degree = degree;
        end
        %%%%%%%%%%%%%%%%%%%
        % fit the data
        function [] = fit(obj,X,Y)
            if nargin == 1
                X = obj.Xdata;
                Y = obj.Xdata;
            else
                obj.Xdata = X;
                obj.Ydata = Y;
            end
            obj.parameters = spap2(obj.knots,obj.degree,X,Y);
            obj.parameters = spap2(newknt(obj.parameters),obj.degree,X,Y);
        end
        
        %%%%%%%%%%%%%%%%%%%
        % eval the logistic function
        function [Y] = fnval(obj,X)
            Y = fnval(obj.parameters,X);
        end
    end
    
end