classdef phytoAaffine < phytoAgeo & sampleable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % my affine transformation
    % an affine transformation can be anywhere and anyway
    % it can be displaced and translated as an object
    % without changing its action
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAaffine(varargin)
            % super constructor - 2 rank fibre and 0 rank base
            obj = obj@phytoAgeo(2,0);
            % set default views
            obj.view_props.props.LineStyle = '--';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAaffine';
            % init point(s)
            if nargin == 1
                data = varargin{1};
                tr = size(data,2);
                data = varargin{1};
                for e = 1:tr
                    putTrial(obj,data(:,e));
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
                op = @(x)phytoAaffine.pullBack(x);
                distrib(obj,fibreI,trialI,op);
            elseif degree > 0
                E = eye(dim(obj)+1);
                op = @(x)phytoAaffine.pushForward(E,x);
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
        % put Trial
        function [] = putTrial(obj,T,trialI,fibreI,toNormalize)
            % set default trial and fibre index
            if nargin == 2;trialI = nTrials(obj)+1;fibreI=1;end
            % set default fibre index
            if nargin == 3;fibreI=1;end
            %%%%%%%%%%%%
            % move to rank-2 object
            T = phytoGeo.toARank2(T);
            %%%%%%%%%%%%
            % get pull back
            nT = phytoAgeo.uNormalize(T);
            %%%%%%%%%%%%
            % object d
            obj.d(trialI,fibreI,:,:) = nT;
            % store pull back
            %obj.d = cat(trialDim,obj.d,shiftdim(nT,insertLevel));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get Trial
        function [t] = getTrial(obj,trialI,fibreI,toShape)
            if nargin == 1;fibreI = 1;end
            % hardcoded problem
            t = obj.d(trialI,fibreI,:,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get subset
        function [s] = getSubSet(obj,trialI,fibreI,toNormalize)
            % set default trial and fibre index
            if nargin == 2;toNormalize = 1;end
            % set default trial and fibre index
            if nargin == 3;fibreI = 1;toNormalize = 1;end
             % set default trial and fibre index
            if nargin == 4;fibreI = 1;toNormalize = 1;end
            % unnormalize the requested subset
            % normalize(obj,1,trialI,fibreI);
            % get trial
            t = getTrial(obj,trialI,fibreI,1);
            % renormalize the requested subset
            if toNormalize;normalize(obj,-1,trialI,fibreI);end
            % cast for return
            s = phytoApoint(t');
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
        % displace
        function [] = displace(obj,v)
            obj.d(1:end-1,end) = obj.d(1:end-1,end) + v.d(1:end-1,end);
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



 %{
        
        %}
%}