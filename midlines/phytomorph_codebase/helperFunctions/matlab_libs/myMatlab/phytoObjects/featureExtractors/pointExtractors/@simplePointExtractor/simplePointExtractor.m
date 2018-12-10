classdef simplePointExtractor  < pointExtractor
    
    properties
        combineFunc;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = simplePointExtractor(varargin)
                % set point extract function
            obj.funcList{1} = @(fM,para)pointsFromMap(fM,para);
                % for non max supports
            obj.para{1}.nonmaxsuppts.rad.value = 35;
            obj.para{1}.nonmaxsuppts.rad.notes = 'radius for finding the local max';
                % for threshold
            obj.para{1}.binaryOperator.op.value = @gt;
            obj.para{1}.binaryOperator.threshold.value = .5;
                % for mask
            obj.para{1}.MSK = [];
                % for combineFunc
            obj.combineFunc = [];
                % override defaults
            if nargin >= 1;  obj.para{1}.nonmaxsuppts.rad.value = varargin{1};end;
            if nargin >= 2;  obj.para{1}.binaryOperator.threshold.value = varargin{2};end;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the mask
        function [MSK] = getMask(obj)
            MSK = obj.para{1}.MSK;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the mask
        function [] = setMask(obj,MSK)
            obj.para{1}.MSK = MSK; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % install combination function
        function [] = setCfunc(obj,cFunc)
            obj.combineFunc = cFunc;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get installed combine func
        function [r] = getCfunc(obj)
            r = obj.combineFunc;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clear combine func
        function [] = clearCfunc(obj)
            obj.combineFunc = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if installed, operate with combine func and follow with
        % extraction
        function pS = extractPoints(obj,X)
            % if combine is given,then run
            if ~isempty(obj.combineFunc)
                X = obj.combineFunc(X);
            end
            % run given function for extraction
            pL = obj.funcList{1}(X,obj.para{1});
            % if isempty
            if isempty(pL);
                pL = [(size(X)+1)/2];
            end
            pS = myHS_X('phytoApoint');
            % loop over points
            for e = 1:size(pL,1)
                % construct set of phytoApoint                        
                pS{e} = phytoApoint([pL(e,:) 1]);            
            end
        end
    end
end