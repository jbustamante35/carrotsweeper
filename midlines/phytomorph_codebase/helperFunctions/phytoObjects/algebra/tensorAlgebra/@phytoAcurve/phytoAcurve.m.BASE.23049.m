classdef phytoAcurve < myTb
    
    properties
        
    end
    
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAcurve(varargin)
            % super constructor - 1 rank fibre and 0 rank base
            obj = obj@myTb();
            obj.setAllRank(1,1);
            % set default views
            obj.view_props.props.LineStyle = '-';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAcurve';
            % init point(s)
            if nargin == 1               
               obj.setData(varargin{1});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,frame,vProps)
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
            
            % reader
            I = myReader(image);               
            affine = rep(obj,[],'phytoAffine');
            affine.normalize(1);
            d = phytoGeo.affineX(domain.d,affine.d);
            % interpolation
            s = myInterp(I,d(1:2,:)');
            % reshape to domain
            s = reshape(s,[domain.gen_parameters.sz size(s,2)]);
            % return patch
            s = phytoPatch(s);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % diff
        function [dc] = diff(obj)
                data = diff(obj,d,1,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate pb_pca
        function [ac] = generate_pbpca(obj,imageStack,Domain)
            for e = 1:size(obj.d,1)
                pt = obj.d(e,:);
                ele = pbPCA_forCurve(pt,imageStack,Domain);
                
                frameSet(e,:,:) = ele.d;
            end
            ac = phytoAsequence(frameSet);
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