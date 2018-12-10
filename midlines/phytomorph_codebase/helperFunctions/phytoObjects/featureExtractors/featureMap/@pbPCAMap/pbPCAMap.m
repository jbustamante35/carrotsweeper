classdef pbPCAMap  < phytoFunc
    
    properties
        
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = cornerFeatureMap(varargin)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init func
                %%%%%%%%%%%%%%%%%%%%%%%%%%%    
            obj.func = @(image,para)cornerMap(image,para);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init notes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.notes = 'function call for generating corner feature map';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init default values for para
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.para.vars.sig.value = 11;
            obj.para.vars.sig.notes = 'sigma for filtering image';
            obj.para.vars.gradPara.sz.value = 0;
            obj.para.vars.gradPara.sz.notes = 'parameters for generating gradient';
            obj.para.vars.gradPara.method.value = 'finite';
            obj.para.vars.gradPara.method.notes = 'differential method';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init non-default parameters
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin >= 1;  obj.para.vars.sig.value = varargin{1};end;            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute feature map
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fM] = computeFeatureMap(obj,X)
            fM = obj.func(X,obj.para);
        end
    end
end