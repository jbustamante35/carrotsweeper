classdef phytoApoint < phytoAgeo & sampleable
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
            % super constructor - 2 rank fibre and 0 rank base
            obj = obj@phytoAgeo(2,0);
            % set default views
            obj.view_props.props.LineStyle = '--';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoApoint';
            % init point(s)
            if nargin == 1
                data = varargin{1};
                tr = size(data,2);
                for e = 1:tr
                    putTrial(obj,data(:,e));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call normalize
        function [] = normalize(obj,degree,trialI,fibreI)
            % if degeee < 0 then create affine transform   
            % if degree > 0 then project into "grand" space            
            % default trials and fibre index
            if (nargin == 2);trialI = 1:obj.nTrials();fibreI = 1;end
            % default fibre index
            if (nargin == 3);fibreI = 1;end;
            % if neg then pull back
            if degree < 0
                op = @(x)phytoApoint.pullBack(x);
                distrib(obj,fibreI,trialI,op);
            elseif degree > 0
                E = eye(dim(obj)+1);
                op = @(x)phytoApoint.pushForward(E,x);
                distrib(obj,fibreI,trialI,op);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation        
        function [r] = rep(obj,frame,type)
            % set frame
            %if isempty(frame);frame = obj.bf;end
            % create eye
            d = eye(size(obj.d,1));
            % switch on rep type
            switch type
                case 'phytoPoint'
                    r = phytoPoint();
                    for e = 1:nTrials(obj)
                        t = getTrial(obj,e);
                        t = phytoAgeo.uNormalize(d,t);
                        t = t(:,end-1);
                        r.putTrial(t);
                    end
                case 'phytoApoint'
                    r = obj.copy();                    
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,frame,vProps,trialI,fibreI)
            % if degeee < 0 then create affine transform   
            % if degree > 0 then project into "grand" space
            % n is trial index : is no n then normalize all
            if (nargin == 3);vProps = [];trialI =  1:obj.nTrials();fibreI = 1;end
            if (nargin == 4);trialI =  1:obj.nTrials();fibreI = 1;end
            if (nargin == 5);fibreI = 1;end;
            % project via frame
            %normalize(obj,1,trialI,fibreI);
            % super view call
            view@phytoAgeo(obj,h,frame,vProps,trialI,fibreI);
            % project via frame
            %normalize(obj,-1,trialI,fibreI);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put trial
        function [] = putTrial(obj,T,trialI,fibreI,toNormalize)
            %%%%%%%%%%%%%%%%%%%%%%%%
            % set default trial and fibre index
            if nargin == 2;trialI = nTrials(obj)+1;fibreI=1;toNormalize = 1;end
            % set default fibre index
            if nargin == 3;fibreI=1;toNormalize = 1;end
            % set toNormalize 
            if nargin == 4;toNormalize = 1;end
            %%%%%%%%%%%%%%%%%%%%%%%%
            % move to rank-2 object
            T = phytoAgeo.toARank2(T);
            %%%%%%%%%%%%
            % get pull back
            if toNormalize;T = phytoAgeo.uNormalize(T);end
            %%%%%%%%%%%%
            % store the data point
            obj.d(trialI,fibreI,:,:) = T;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get Trial
        function [t] = getTrial(obj,trialI,fibreI,toShape)
            % set default trial and fibre index
            if nargin == 2;fibreI = 1;toShape = 0;end
            % set default fibre index
            if nargin == 3;toShape = 0;end
            % hardcoded problem
            t = obj.d(trialI,fibreI,:,:);
            % squeeze
            %if ~toShape;t = squeeze(t);end
            if ~toShape;t = shiftdim(t,2);end
            % toShape
            if toShape;t = t(:,:,:,end);sz = size(t);t = reshape(t,[sz(1),sz(end)]);t = t';end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get subset
        function [s] = getSubSet(obj,trialI,fibreI,toNormalize)
            % set default trial and fibre index
            if nargin == 2;fibreI = 1;toNormalize = 1;end
            % set default trial and fibre index
            if nargin == 3;toNormalize = 1;end            
            % unnormalize the requested subset
            % normalize(obj,1,trialI,fibreI);
            % get trial
            t = getTrial(obj,trialI,fibreI,1);
            % renormalize the requested subset
            if toNormalize;normalize(obj,-1,trialI,fibreI);end
            % cast for return
            typeC = str2func(class(obj));
            s = typeC(t);
            %s = phytoApoint(t');
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % distribute operator over the fibre and trial vectors
        function [] = distrib(obj,fibreI,trialI,op)
            % over each fibreI and each trialI
            for f = 1:numel(fibreI)
                for t = 1:numel(trialI)
                    tmp = getTrial(obj,trialI(t),fibreI(f));
                    tmp = op(tmp);
                    putTrial(obj,tmp,trialI(t),fibreI(f),0);
                end
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dimension of vector or point
        function [r] = dim(obj)
            r = size(obj.d,4)-1
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
    end
end