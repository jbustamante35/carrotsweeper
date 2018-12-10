classdef continuousDomainMap
    % a set of domains and a map - see cordinateTransform
    properties
        domainArray;
        mapFunc;
    end
    
    methods
        function [obj] = continuousDomainMap(mapFunc,domainArray)
            obj.domainArray = domainArray;
            obj.mapFunc = mapFunc;
        end
        
        function [] = setDomain(obj,domain,n)
            obj.domainArray{n} = domain;
        end
        
        function [domain] = getDomain(obj,n)
            domain = obj.domainArray{n};
        end
        
        function [] = deleteDomain(obj,n)
            obj.domainArray(n) = [];
        end
        
        function [] = insetDomain(obj,domain,n)
            obj.domainArray{n} = domain;
        end
        
        function [domain] = evaluateAt(obj,para,n)
            if nargin == 3
                domain{1} = obj.mapFunc(obj.domainArray{n},para);
            else
                for e = 1:numel(obj.domainArray)
                    domain{e} = obj.mapFunc(obj.domainArray{e},para);
                end
            end
        end
        
        function [] = plot(obj,skip,CL,P)
            for e = 1:numel(obj.domainArray)
                obj.mapFunc.plot(skip,CL,obj.domainArray{e},P);
            end
        end
        
    end
    
end

%{

    [X1(:,:,1),X1(:,:,2)] = ndgrid(linspace(0,6,10),linspace(-pi,pi,100));
    [X2(:,:,1),X2(:,:,2)] = ndgrid(linspace(6,12,10),linspace(-pi,pi,100));
    [X3(:,:,1),X3(:,:,2)] = ndgrid(linspace(12,24,10),linspace(-pi,pi,100));
    domainArray{1} = X1;
    domainArray{2} = X2;
    domainArray{3} = X3;

    func = @(X,P)cat(3,P(2)*cos(P(1))*X(:,:,1).*cos(X(:,:,2))-P(3)*sin(P(1))*X(:,:,1).*sin(X(:,:,2))+P(4),P(2)*sin(P(1))*X(:,:,1).*cos(X(:,:,2))+P(3)*cos(P(1))*X(:,:,1).*sin(X(:,:,2))+P(5));
    continuousMapper = cordinateTransform(func);
    para = [0 1 1 100 100];

    domainMapper = continuousDomainMap(continuousMapper,domainArray);
    d = domainMapper.evaluateAt(para);

    domainMapper.plot([3 3],{'r'},para)




%}