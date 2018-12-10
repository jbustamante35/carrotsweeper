classdef phytoApatchSequence < myTb
    
    properties
        
    end
    
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoApatchSequence(varargin)
            % super constructor - 1 rank fibre and 0 rank base
            obj = obj@myTb();
            obj.setAllRank(1,2);
            % set default views
            obj.view_props.props.LineStyle = '-';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoApatchSequence';
            % init point(s)
            if nargin == 1               
               obj.setData(varargin{1});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,vProps)
            S.type = '()';
            S.subs = {':'};            
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);            
            tensorView(subsref(obj,S),uProps,h);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation  - not done!!      
        function [r] = rep(obj,frame,type)
            switch type
                case 'phytoApoint'
                   
                case 'phytoAaffine'
                   
                case 'phytoAcurve'
                         
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reparameterize
        function [] = arcLength(obj,type,value)
            obj.d = arcLength(obj.d,type,value);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % length
        function [l] = length(obj)
            l = length(obj.d)';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % length
        function [k] = kurvature(obj,smoothValue)
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pull back
        function [d] = pullBack(d)
            for e = 1:size(d,1)
                d(e,:,:) = inv(d(e,:,:));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push forward
        function [d] = pushForward(d,v)
            for e = 1:size(d,1)
                d(e,:,:) = d(e,:,:)*v;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct from A point
        function [cur] = contructFromApoint(p)
            d = p();
            d = [d(1:end-1) 1 d(end)];
            cur = phytoAcurve(d);
        end
    end
end