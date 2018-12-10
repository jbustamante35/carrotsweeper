classdef continuousMap < handle
    
    properties
        
        
        multiLinMap;
        mapper;
        
        
        domainFunc;
    end
    
    methods
        function [obj] = continuousMap(multiLinMap,mapper,domainFunc)
            obj.multiLinMap = multiLinMap;
            obj.mapper = mapper;
            obj.domainFunc = domainFunc;
        end
          
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'{}')
                try 
                    fprintf(['error in index-{} for continuousMap map.\n']);
                catch ME
                    fprintf(['error in index-{} for continuousMap map.\n']);
                end
            elseif strcmp(S(1).type,'.')
                try
                    fprintf(['error in index-. for continuousMap map.\n']);
                catch
                    fprintf(['error in index-. for continuousMap map.\n']);
                end
            elseif strcmp(S(1).type,'()')
                try
                    % get the domains at the para
                    domain = obj.domainFunc.evaluateAt(S(1).subs{:});
                    % for each domain
                    for e = 1:numel(domain)
                        varargout{1}{e} = obj.mapper(obj.multiLinMap.M,domain{e}(:,:,1),domain{e}(:,:,2));
                    end
                catch ME
                    fprintf(['error in index-() for continuousMap map.\n']);
                end
            else ME
                fprintf(['error in index-ANY for multilinear map.\n']);
            end
        end
        
        
        
    
    end
end
%{
    I = imread('/home/nate/Downloads/20170329n02_04.tif');
    mI = multiLinearMap(I);

    

    
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


    M = continuousMap(mI,str2func('ba_interp2'),domainMapper);
    sample = M(para);
    imshow(I,[]);
    hold on
    cM.plot([3 3],{'r'},Xu,para);
    imshow(sample,[]);












%}