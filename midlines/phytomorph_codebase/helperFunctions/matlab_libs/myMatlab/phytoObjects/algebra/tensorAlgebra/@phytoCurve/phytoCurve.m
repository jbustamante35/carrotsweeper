classdef phytoCurve < phytoAgeo
    
    properties
        
    end
    
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoCurve(varargin)
            % super constructor
            obj = obj@phytoAgeo(2,1);            
            % set default views
            obj.view_props.props.LineStyle = '-';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoCurve';
            % init point(s)
            if nargin == 1
                data = varargin{1};
                tr = size(data,3);
                for e = 1:tr
                    putTrial(obj,data(:,:,e));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert trial
        function [] = putTrial(obj,T)
            T = phytoCurve.toAffine(T);
            putTrial@phytoAcurve(obj,T);
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % toAffine
    methods (Static)
        function [data] = toAffine(data)
            data = cat(1,data,ones(1,size(data,2),size(data,3)));
        end
    end

end

%{

    %%%%%%%%%%%%%%%%
    % graph
    G = myG();

    %%%%%%%%%%%%%%%%
    % generate X rand nodes
    N = 100;
    for e = 1:N
        tN = myN();
        pt = phytoPoint(rand(2,1));
        pt.normalize();
        G.putNode(tN);
    end

    ten = myT();
%}