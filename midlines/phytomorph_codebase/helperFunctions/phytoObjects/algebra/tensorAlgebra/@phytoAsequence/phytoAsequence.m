classdef phytoAsequence < myTb
    
    properties
        
    end
    
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAsequence(varargin)
            % super constructor - 2 rank fibre and 1 rank base
            obj = obj@myTb();
            obj.setAllRank(1,2);
            % set default views
            obj.view_props.props.LineStyle = '-';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAsequence';
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample image at domain
        function [s] = sample(obj,image,domain)
            try
                % sample over the affine objects at the domain
                for e = 1:baseSize(obj)
                    S.type = '()';
                    S.subs{1} = e;
                    ele = squeeze(subsref(obj,S));
                    affine = phytoAaffine(ele);
                    ele = affine.sample(image,domain);
                    data(e,:,:) = ele.d;
                end                
            catch ME
               ME; 
            end
            s = phytoApatchSequence(data);
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % angle
        function [data] = angle(obj)
            % sample over the affine objects at the domain
            for e = 1:baseSize(obj)
                S.type = '()';
                S.subs{1} = e;
                ele = squeeze(subsref(obj,S));
                affine = phytoAaffine(ele);                
                ele = affine.angle();
                data(e) = ele;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % velocity
        function [data] = velocity(obj)
            path = [];
            % sample over the affine objects at the domain
            for e = 1:baseSize(obj)
                S.type = '()';
                S.subs{1} = e;
                ele = squeeze(subsref(obj,S));
                path = [path ele(:,end)];
            end
            dP = diff(path,1,2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % velocity
        function [] = arcLengthPara(obj,baseCurve)            
            paraL = baseCurve.iLength();
            N = round(paraL(end));
            L = linspace(paraL(1),paraL(end),N);            
            obj.d = interp1(paraL,obj.d,L);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % path
        function [path] = path(obj)
            path = [];
            % sample over the affine objects at the domain
            for e = 1:baseSize(obj)
                S.type = '()';
                S.subs{1} = e;
                ele = squeeze(subsref(obj,S));
                path = [path ele(:,end)];
            end
            path = phytoAcurve(path(1:2,:)');
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