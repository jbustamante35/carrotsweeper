classdef phytoAgeo < myT & geometrical
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct for geo-objects : 
    % [base_trial,fibre_trial,base_manifold,fibre_manifold]
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        boundingBox;    % boundingbox
        rank;           % obj rank [fibre,base] rank - rank for base and fibre - fixed at 2 for fibre
        sz;             % size of data
        di;             % dictionary object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        base_trial_key = 'tr_b';
        fibre_trial_key = 'tr_f';
        base_manifold_key = 'mi_b';
        fibre_manifold_key = 'mi_f';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAgeo(varargin)
            obj@geometrical();
            obj@myT();
            %%%%%%%%%%%%%%%%
            % bounding box
            obj.boundingBox = [];
            obj.rank = 0;
            obj.di = dictionary();
            %%%%%%%%%%%%%%%%
            % trial base
            obj.di.insertTerm(obj.base_trial_key);
            % trial fibre
            obj.di.insertTerm(obj.fibre_trial_key);
            %%%%%%%%%%%%%%%%
            % init dictionary
            if nargin == 2
                % set geo-rank
                obj.rank = [varargin{1} varargin{2}];
                % fill manifold
                for r = 1:obj.rank(2)
                    obj.di.insertTerm([obj.base_manifold_key num2str(r)]);
                end
                 % fill extra d
                for r = 1:(obj.rank(1))
                    obj.di.insertTerm([obj.fibre_manifold_key num2str(r)]);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call normalize
        function [] = normalize(obj,degree,trialI,fibreI)
            % if degeee < 0 then create affine transform   
            % if degree > 0 then project into "grand" space            
            % default trials to all and fibre index to base fibre
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
        % number of base spaces
        function [n] = nTrials(obj)
            r = 1;
            n = size(obj.d,r);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of fibre spaces
        function [n] = nFibres(obj)
            r = 2;
            n = size(obj.d,r);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get trial dim
        function [n] = baseTdim(obj)
            n = obj.di.lookUp(obj.base_trial_key);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get fibre trial dim
        function [n] = fibreTdim(obj)
            n = obj.di.lookUp(obj.fibre_trial_key);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get rank of fibre
        function [n] = fibreRank(obj)
            n = obj.rank(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % base rank  
        function [n] = baseRank(obj)
            n = obj.rank(2);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % base dims
        function [r] = baseDims(obj)
            r = obj.di.seachTerm(obj.base_manifold_key);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fibre dims
        function [r] = fibreDims(obj)
            sz = size(obj.d);
            r = sz(end-1:end);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to affine stream
        function [] = toFibreStream(obj,formRank)
            % set default stream rank
            if (nargin == 1)
                formRank = 2;
            end
            % store size of data
            obj.sz = size(obj.d);
            % get size of data
            sz = size(obj.d);
            % get the base
            r = baseDims(obj);
            % prod of base dims
            nD = prod(sz(r));
            % create return
            nSz = sz;
            % fill in size
            if ~isempty(r)
                nSz = [nSz(1:(r(1)-1)) nD nSz(r(end)+1:end)];
            end
            % output rank
            if (formRank == 1)
                nSz = [nSz prod(nSz(end-1:end))];
            end
            % perform reshape of object
            dReshape(obj,nSz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to orginal shape
        function [] = toOshape(obj)
            dReshape(obj,obj.sz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dualReshape
        function [] = dReshape(obj,sz)
            obj.d = reshape(obj.d,sz);
            obj.bf = reshape(obj.bf,sz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [obj] = subsasgn(obj,S,B)
            % call assign via super
            subsasgn@myT(obj,S,B);
            % if normalize flag is true
            if obj.s;obj.normalize(-1);end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,frame,vProps,trialI,fibreI)
            % if not defined, then set frame to default
            % if isempty(frame);frame=obj.bf;end
            % init vProps is non def
            if nargin == 3;vProps = [];end
            % get subset
            s = getSubSet(obj,trialI,fibreI);
            % get subset
            % subS = getTrial(obj,trialI,fibreI);
            % project via frame         
            s.normalize(1);
            % call view with projected data
            view@myT(s,h,frame,vProps);
            % project back            
            % normalize(obj,-1,trialI,fibreI);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate bounding box
        function [] = generateBoundingBox(obj)
            mM = min(obj.d,2);
            MM = max(obj.d,2);
        end
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % xform via affine verb
        function [d] = affineX(d,verb)
            if size(verb) == (size(d,1) + 1)
                d = [d;ones(1,size(d,2))];
            end
            d = verb*d;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % "universal" normalization
        function [r] = uNormalize(data,verb)
            % expected form is stream of 2nd rank tensors
            r = zeros(size(data));
            if nargin == 1
                % iterate over elements
                for e = 1:size(data,3)
                    r(:,:,e) = inv(data(:,:,e));
                end
            else
                % iterate over elements
                for e = 1:size(data,3)
                    r(:,:,e) = verb(:,:,e)*data(:,:,e);
                end 
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % toRank2
        function [r] = toARank2(data)
            S = size(data);
            if (ndims(data) == 1) | (S(1) ~= S(2))
                r = eye(size(data,1));
                r(:,end) = data;
            elseif ndims(data) == 2
                r = data;
            end
        end
    end
    
    methods (Abstract)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % implement normalize
        %[] = normalize(obj,degree,trialI,fibreI);
    end
end