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
                if ~isjava(varargin{1});
                    obj.setData(varargin{1});
                else
                    obj = phytoAcurve.fromBson(varargin{1});
                end
            end
        end
        % toBson
        function [xferO] = toBson(obj,xferO)
            import phytoG.locked.BdataObjects.geometry.implementations.*;    
            if (nargin == 1);xferO = phytoAcurve();end
            xferO = toBson@myTb(obj,xferO);
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
        % representation - NOT DONE!  
        function [r] = rep(obj,frame,type)
            switch type
                case 'phytoApoint'
                   
                case 'phytoAaffine'
                   
                case 'phytoAcurve'
                         
            end
        end
        % sample along curve
        function [s] = sample(obj,image,domain)
            for e = 1:size(obj.d,1)
                affineData = obj(e);
                phytoAaffine(affineData);
                sample(obj,image,domain) 
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % diff
        function [dc] = diff(obj)
                data = diff(obj,d,1,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate pb_pca
        function [ac] = generate_pbpca(obj,imageStack,Domain,extraDims)
            if nargin == 3;extraDims = [];end
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ilength
        function [L] = iLength(obj)
            dL = diff(obj.d,1,1);
            L = sum(dL(:,1:2).*dL(:,1:2),2).^.5;
            L = cumsum([0;L]);
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % acrlength parameterize
        function [k] = arcLengthPara(obj)
            obj.d = arcLength(obj.d,'arcLen',[]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % tangent space
        function [TM] = generateTangentSpace(obj,level)  
            TM = igetFrame(obj.d(:,1:2),level);
            TM = cat(3,TM,obj.d(:,1:2));
            zTM = zeros(size(TM,1),3,3);
            for e = 1:size(TM,1)
                zTM(e,:,:) = eye(3);
                zTM(e,1:2,:) = TM(e,:,:);
            end
            TM = phytoAsequence(zTM);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % tangent space @ location
        function [TM] = generateTangentSpace_atP(obj,point,level)
            TM = igetFrame_atP(obj.d(:,1:2),point,level);            
            zTM = eye(3);
            zTM(1:2,1:2) = TM;
            zTM = shiftdim(zTM,-1);
            TM = phytoAsequence(zTM);
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
        
        % fromBson
        function [out] = fromBson(varargin)
            in = varargin{1};
            
            if (nargin == 1);out = phytoAcurve();
            else out = varargin{2};end
            
            out = myTb.fromBson(in,out);
        end
    end
end
