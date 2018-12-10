classdef phytoAcurve < phytoAgeo
    
    properties
        
    end
    
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAcurve(varargin)
            % super constructor - fibreRank=2,baseRank=1;
            obj = obj@phytoAgeo(2,1);            
            % set default views
            obj.view_props.props.LineStyle = '-';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAcurve';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation  - not done!!      
        function [r] = rep(obj,frame,type)
            % set frame
            %if isempty(frame);frame = obj.bf;end            
            % switch on rep type
            switch type
                case 'phytoApoint'
                    r = phytoGeo.affineX(inv(obj.d),frame);
                    r = phytoPoint(mean(r,2));
                case 'phytoAaffine'
                    tpt = phytoGeo.affineX(obj.d,inv(frame));
                    r = phytoAffine(eye(size(obj.d,1)));                    
                    r(1:2,3) = tpt(1:2);
                case 'phytoAcurve'
                    r = phytoGeo.affineX(inv(obj.d),frame);
                    r = phytoPoint(mean(r,2));        
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
             % set default trial and fibre index
            if nargin == 2;trialI = nTrials(obj)+1;fibreI=1;end
            % set default fibre index
            if nargin == 3;fibreI=1;end
            nP = size(T,2);
            objDim = size(T,1);
            %%%%%%%%%%%%
            % normalize single curve
            curveD = zeros(nP,objDim,objDim);
            for cP = 1:nP
                pt = T(:,cP);
                pt = phytoGeo.toARank2(pt);
                curveD(cP,:,:) = shiftdim(phytoAgeo.uNormalize(pt),-1);
            end
            obj.d(trialI,fibreI,:,:,:) = curveD;
            %{
            %%%%%%%%%%%%
            % shift dim
            curveBF = shiftdim(curveBF,insertLevel);
            curveD = shiftdim(curveD,insertLevel);
            %%%%%%%%%%%%
            % store
            obj.d = cat(trialDim,obj.d,curveD);
            obj.bf = cat(trialDim,obj.bf,curveBF);
            %}
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get trial
        function [t] = getTrial(obj,trialI,fibreI,toShape)
            % set defaults for trialI,fibreI and toShape
            if nargin == 1;trialI = 1;fibreI = 1;toShape = 0;end
            % set default trial and fibre index
            if nargin == 2;fibreI = 1;toShape = 0;end
            % set default fibre index
            if nargin == 3;toShape = 0;end
            % hardcoded problem third from end is along curve
            t = obj.d(trialI,fibreI,:,:,:);
            % shift to release the fibre and trial dim
            if ~toShape;t = shiftdim(t,2);end
            % toShape
            if toShape;t = t(:,:,:,end);sz = size(t);t = reshape(t,[sz(1),sz(end)]);end
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
            s = phytoAcurve(t');
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
        % displace
        function [] = displace(obj,v)
            obj.d(1:end-1,end) = obj.d(1:end-1,end) + v.d(1:end-1,end);
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