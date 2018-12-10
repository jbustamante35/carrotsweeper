classdef phytoGeo < phytoAgeo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct for geo-objects
    % todo: remove the trial count
    properties
       
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoGeo(varargin)
            obj@phytoAgeo();
        end
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call normalize
        function [] = normalize(obj,h,frame,vProps,trialI,fibreI)
            % if degeee < 0 then create affine transform   
            % if degree > 0 then project into "grand" space            
            if (nargin == 2);trialI = 1:obj.nTrials();fibreI = 1;end
            if (nargin == 3);fibreI = 1;end;
            % if neg then pull back
            if degree < 0
                op = @(x)phytoAcurve.pullBack(x);
                distrib(obj,fibreI,trialI,op);
            elseif degree > 0
                E = eye(dim(obj)+1);
                op = @(x)phytoAcurve.pushForward(E,x);
                distrib(obj,fibreI,trialI,op);
            end
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % toRank2
        function [r] = toARank2(data)      
            r = eye(size(data,1));
            r(:,end) = data;
        end
    end
    
    methods (Abstract)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % implement normalize
        [] = normalize(obj,degree);
    end
end