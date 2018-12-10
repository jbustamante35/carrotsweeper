classdef gsS < myHS
    
    % phyto-graph-sample-set - gss
    properties
        E;  % edges
        N;  % nodes
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = gsS(varargin)
            % call super constructor
            obj = obj@myHS();
            % init two constraint sets which contain patches
            obj.E = myHS_X('phytoPatch');
            obj.N = myHS_X('phytoPatch');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert phytoPatch
        function [obj] = insertSample(obj,patch,setName)
            obj.(setName).putElement(patch);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert phytoPatch - edge
        function [obj] = insertEdge(obj,patch,setName)
            insertSample(obj,patch,'E');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert phytoPatch - node
        function [obj] = insertNode(obj,patch,setName)
            insertSample(obj,patch,'N');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert sample
        function [obj] = view(obj,h,frame,varargin)
            
            for e = 1:numel(obj.E)
                p = obj.E{e};
                p.view(h,frame);
                drawnow;                
                waitforbuttonpress;
            end
            
            for e = 1:numel(obj.N)
                p = obj.N{e};
                p.view(h,frame);
                drawnow;
                waitforbuttonpress;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation
        function [r] = rep(obj,frame,type)
        
        end
        
    end
    
    
    
    
end