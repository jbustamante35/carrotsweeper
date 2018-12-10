classdef phytoDomain < phytoGeo
    
    properties
        gen_parameters;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoDomain(varargin)
            % super constructor
            obj = obj@phytoGeo();
            % set as defined object
            obj.setState(0);
            % set view props
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAcurve';
            
            if nargin == 1
                obj.gen_parameters = varargin{1};
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate
        function [] = generateDomain(obj)
            if isempty(obj.d)
                % call to generate domain
                obj.gen_parameters = createDomain(obj.gen_parameters);
                % store domain data
                obj.d = obj.gen_parameters.d;
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
        % call normalize
        function [] = normalize(obj,degree)
            % if degee < 0 then create affine transform            
            if degree < 0
                obj.bf = [eye(3)];
                obj.bf(1:2,3) = -mean(obj.d,2);
                obj.d = phytoGeo.affineX(obj.d,obj.bf);
            elseif degree > 0
                if ~isempty(obj.bf)
                    obj.d = phytoGeo.affineX(inv(obj.d),obj.bf);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,vProps)           
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);
            % call view with projected data
            view@myT(obj,h,vProps);
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