classdef phytoPoint < phytoFB%< phytoApoint
    
    properties
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoPoint(varargin)
            % super constructor
            obj = obj@phytoFB();
            % set default views
            obj.view_props.props.LineStyle = '--';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoPoint';
            % init point(s)
            if nargin == 1
                tr = size(varargin{1},2);
                data = varargin{1};
                for e = 1:tr
                    putTrial(obj,data(:,e));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert trial
        function [] = putTrial(obj,T,trialI,fibreI,toNormalize)
            %%%%%%%%%%%%%%%%%%%%%%%%
            % set default trial and fibre index
            if nargin == 2;trialI = nTrials(obj)+1;fibreI=1;toNormalize = 1;end
            % set default fibre index
            if nargin == 3;fibreI=1;toNormalize = 1;end
            % set toNormalize 
            if nargin == 4;toNormalize = 1;end
            %%%%%%%%%%%%%%%%%%%%%%%%
            % normalize here means pullBack and affine
            if toNormalize;T = phytoPoint.toAffine(T);end
            %%%%%%%%%%%%%%%%%%%%%%%%
            % put trial
            putTrial@phytoApoint(obj,T,trialI,fibreI,toNormalize);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get trial
        function [t] = getTrial(obj,trialI,fibreI,toShape)
            % set default trial and fibre index
            if nargin == 2;fibreI = 1;toShape = 0;end
            % set default fibre index
            if nargin == 3;toShape = 0;end
            % get trial information
            t = getTrial@phytoApoint(obj,trialI,fibreI,toShape);            
            % toShape
            if toShape;t = t(1:end-1,:);end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % toAffine
    methods (Static)
        function [data] = toAffine(data)
            data = [data;ones(1,size(data,2))];
        end
    end
end
