classdef tracer < lManifold
    
    properties
        xSZ;
        ySZ;
    end
    
    methods
        function [obj] = tracer()
            obj@lManifold();
        end
        
        function [] = setXsz(obj,sz)
            obj.xSZ = sz;
        end
        
        function [] = setYsz(obj,sz)
            obj.ySZ = sz;
        end
        
        function [] = viewDomain(obj)            
            for e = 1:size(obj.rawX,1)
                patch = reshape(obj.rawX(e,:),obj.xSZ);
                imshow(patch,[])
                drawnow
            end
        end
        
        function [] = viewManifoldBase(obj)
            M = [];
            for e = 1:obj.nGroups
                M = cat(4,M,reshape(obj.Ux(e,:),obj.xSZ));
            end
            montage(M)
        end
        
        
        function [] = viewCoandDo(obj)
            for e = 1:size(obj.rawX,1)
                patch = reshape(obj.rawX(e,:),obj.xSZ);
                hold off
                imshow(patch,[]);   
                hold on
                curve = reshape(obj.rawY(e,:),obj.ySZ);
                curve = bsxfun(@plus,curve,(size(patch)-1)/2);
                plot(curve(:,1),curve(:,2),'r');
                drawnow
            end
        end
        
    end
end