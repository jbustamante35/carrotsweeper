classdef phytoAdomain < myTb %< phytoGeo
    
    properties
        gen_parameters;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAdomain(varargin)
            % super constructor
            obj = obj@myTb();
            % set base and fibre rank
            obj.setAllRank(0,2);
            % set view props
            obj.view_props.props.LineStyle = '-';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAcurve';
            % set view props
            if nargin == 1
                obj.gen_parameters = varargin{1};
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,vProps)
             % view affine
            S.type = '()';
            S.subs = {}; 
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);            
            tensorView(obj.view_rep,uProps,h);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate
        function [] = generateDomain(obj)
            if isempty(obj.d)
                % call to generate domain
                obj.gen_parameters = createDomain(obj.gen_parameters);
                % set data
                setData(obj,obj.gen_parameters.d);
                % remove from 
                obj.gen_parameters.d = [];
                % store representation data
                obj.view_rep = obj.gen_parameters.rep';
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate domain
        function [] = clearDomain(obj)
            obj.d = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation        
        function [r] = rep(obj,frame,type)
            % set frame
            if isempty(frame);frame = obj.bf;end            
            % switch on rep type
            switch type
                case 'phytoPoint'
                    r = phytoGeo.affineX(obj.d,inv(frame));
                    r = phytoPoint(r(1:2,3));
                case 'phytoAffine'
                    tpt = phytoGeo.affineX(obj.d,inv(frame));
                    r = phytoAffine(eye(size(obj.d,1)));                    
                    r(1:2,3) = tpt(1:2);
                case 'phytoCurve'
                    r = phytoGeo.affineX(inv(obj.d),frame);
                    r = phytoPoint(mean(r,2));        
            end
        end
    end
end

%{
    para.type = 'disk';
    para.value{1} = [0 100 100];
    para.value{2} = [-pi pi 200];
    
    domain = phytoDomain(para);
    domain.generateDomain();
    h = figure;
    domain.view(h,[]);

%}