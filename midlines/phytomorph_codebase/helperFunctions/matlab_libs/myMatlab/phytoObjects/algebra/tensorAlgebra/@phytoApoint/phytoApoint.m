classdef phytoApoint < myTb & persistable %< phytoFB %< phytoAgeo & sampleable
    % my affine transformation
    % an affine transformation can be anywhere and anyway
    % it can be displaced and translated as an object
    % without changing its action
    properties
        
    end
    
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoApoint(varargin)           
            % super constructor - 1 rank fibre and 0 rank base
            obj = obj@myTb();
            obj.setAllRank(0,1);
            % set default views
            obj.view_props.props.LineStyle = '--';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoApoint';
            % init point(s)
            if nargin == 1
                if ~isjava(varargin{1});
                    obj.setData(varargin{1});
                else
                    obj = phytoApoint.fromBson(varargin{1});
                end
            end
        end
        % toBson
        function [xferO] = toBson(obj,xferO)
            import phytoG.locked.BdataObjects.geometry.implementations.*;    
            if (nargin == 1);xferO = phytoApoint();end
            xferO = toBson@myTb(obj,xferO);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,vProps)
            S.type = '()';
            S.subs = {};            
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);            
            tensorView(subsref(obj,S),uProps,h);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation        
        function [r] = rep(obj,frame,type)
            % set frame
            %if isempty(frame);frame = obj.bf;end
            % switch on rep type
            r = [];
            switch type
                case 'phytoApoint'
                    r = obj.copy();
                case 'phytoAaffine'
                    col = obj.d;
                    d = eye(numel(col));
                    d(:,end) = col;
                    r = phytoAaffine(d);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample image at domain
        function [s] = sample(obj,image,domain)
            affine = rep(obj,[],'phytoAaffine');
            % load para for q-reader
            para{1} = affine.d;
            para{2} = domain.d;
            sz = domain.gen_parameters;
            para{3} = sz.sz;
            I = myReader(image,'iatP',para);
            % return patch
            s = phytoApatch(I);
        end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate pb_pca
        function [r] = generate_pbpca(obj,imageStack,Domain,extraDims)
            if nargin == 3;extraDims = [];end
            pt = obj.d;
            r = pbPCA_forCurve(pt,imageStack,Domain,extraDims);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dimension of vector or point
        function [r] = dim(obj)
            r = size(obj.d,4)-1
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure length
        function [l] = norm(obj)
            l = obj.d.*obj.d;
            l = sum(l).^.5;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pull back
        function [id] = pullBack(d)
            id = inv(d);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push forwards 
        function [ret] = pushForward(d,v)
            ret = d*v;
        end
        % fromBson
        function [out] = fromBson(varargin)
            in = varargin{1};
            
            if (nargin == 1);out = phytoApoint();
            else out = varargin{2};end
            
            out = myTb.fromBson(in,out);
        end            
    end
end