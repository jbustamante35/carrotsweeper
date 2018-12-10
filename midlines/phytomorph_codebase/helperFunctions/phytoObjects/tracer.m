classdef tracer < matlab.mixin.Copyable
    properties
        position;
        direction;
        nHood;
        rho;
        rad;
        density;
        
        image;
    end
    
    methods 
        function [obj] = tracer()
        
        end
        
        function [] = setNhoodRho(rho)
            obj.rho = rho;
        end
        
        function [] = setNhoodRad(rad)
            obj.rad = rad;
        end
        
        function [] = setNhoodDensity(density)
            obj.density = density;            
        end
        
        function [] = setPosition(X)
            obj.position = X;
        end
        
        function [] = setDirection(X)
            obj.direction = X;
        end
        
        function [] = setImage(I)
            obj.image = I;
        end
        
        function [] = sampleImage()
            % rotate to last frame
            tmpH = obj.direction(:,:,end)*obj.nHood;
            % displace to last position
            tmpH = bsxfun(@plus,tmpH,obj.position(:,end));
        end
        %{
        function [] = generateH()
            [rho rad] = ndgrid(linspace(0,obj.rho,obj.density(1)),linspace(obj.rad,obj.rad,obj.density(2)));
            obj.nHood = [rho(;).*cos(rad);rho(;).*sin(rad)];
        end
        %}
        
    end
end

%{

T = tracer();
T.setNhoodRho(30);
T.setNhoodRad(pi);
T.setNhoodDensity([30 300]);

%}