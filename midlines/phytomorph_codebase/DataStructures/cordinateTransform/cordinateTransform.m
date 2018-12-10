classdef cordinateTransform
    % maps a ndgrid space to another space with a smooth function
    properties
        map;
    end
    
    methods
        function [obj] = cordinateTransform(map)
            obj.map = map;
        end
        
        function [] = plot(obj,skip,CL,X,P)
            if nargin == 5
                X = obj.map(X,P);
            end
            
            if numel(CL) == 1
                CL{2} = CL{1};
            end
            
            for e1 = 1:skip(1):size(X,2)
                plot(X(:,e1,1),X(:,e1,2),CL{1})
                hold on
            end
            plot(X(:,end,1),X(:,end,2),CL{1})
            
            
            for e1 = 1:skip(2):size(X,1)
                plot(X(e1,:,1),X(e1,:,2),CL{2})
                hold on
            end
            plot(X(end,:,1),X(end,:,2),CL{2})
    
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
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                catch
                    fprintf(['error in index-. for continuousMap map.\n']);
                end
            elseif strcmp(S(1).type,'()')
                try
                    varargout{1} = obj.map(S(1).subs{1},S(1).subs{2});
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


    cT = @(X,P)cat(3,P(2)*cos(P(1))*X(:,:,1).*cos(X(:,:,2))-P(3)*sin(P(1))*X(:,:,1).*sin(X(:,:,2))+P(4),P(2)*sin(P(1))*X(:,:,1).*cos(X(:,:,2))+P(3)*cos(P(1))*X(:,:,1).*sin(X(:,:,2))+P(5));
    


    cM = cordinateTransform(cT);

   

    [X1(:,:,1),X1(:,:,2)] = ndgrid(linspace(0,6,10),linspace(-pi,pi,100));
    [X2(:,:,1),X2(:,:,2)] = ndgrid(linspace(6,12,10),linspace(-pi,pi,100));
    [X3(:,:,1),X3(:,:,2)] = ndgrid(linspace(12,24,10),linspace(-pi,pi,100));


    para = [pi/8 10 5 0 0];
    F = cM(Xo,para);
    

    cM.plot([3 3],{'r'},X1,para);
    cM.plot([3 3],{'g'},X2,para);
    cM.plot([3 3],{'b'},X3,para);



%}