classdef smoothImage  < smoothTensor
    properties
        
    end
    
    methods
        function [obj] = smoothImage(varargin)
            obj = obj@smoothTensor(varargin{:});
        end
        
        function [r] = evaluate(obj,varargin)
            
            varargin{2} = obj;
            
            
            for s = 1:numel(varargin)
                for l = 1:size(varargin{s}.M,3)
                    X{s}{l} = varargin{s}.M(:,:,l);
                end
            end
            
            r = obj.out_class('','',tensor(obj.in_map(X{:})));
        end
        
        function [] = plot(obj)
            imshow(obj.M,[])
            hold on
        end
    end
end






