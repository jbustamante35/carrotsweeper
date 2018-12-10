classdef smoothTensor < tensor
    
    properties
        in_map;
        out_class;
    end
    
    methods
        function [obj] = smoothTensor(in_map,out_class,varargin)
            % if this smooth tensor is attached to any discrete tensor
            obj = obj@tensor(varargin{:});
            obj.in_map = in_map;
            obj.out_class = out_class;
        end
    end
    
    methods
        function [r] = evaluate(obj,varargin)
            
            for s = 1:numel(varargin)
                X{s} = mat2cell(varargin{s}',ones(numel(varargin{s}),1));
                for l = 1:numel(X{s})
                    X{s}{l} = X{s}{l}.ten.M;
                end
            end
            r(e) = obj(e).out_class('','',tensor(obj(e).in_map(X{:})));
            
            
        end
    end
end

%{
 


    











%}