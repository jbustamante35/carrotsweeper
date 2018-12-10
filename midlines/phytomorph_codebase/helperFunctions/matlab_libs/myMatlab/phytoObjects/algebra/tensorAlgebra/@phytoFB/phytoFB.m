classdef phytoFB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct for fibreBundle:
    % during walk home, after concept meeting.
    % i choose this as an implementation. this
    % needs to be as such because the dim of the fibre
    % may not be homogenous over the base. for example
    % the same base of rank 1 such as a curve. will have 
    % the base fibre of the motion of the id over n-dim. 
    % but there could and will be a curve which will have the base fibre
    % of dims 2 (via round up) and a sampling along the fibre of the image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data structure
        gO;             % data structure X(baseT,fibreT).b(baseI).f(fibreI)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bR;             % base rank homogenous over bundle set
        fR;             % set of fibre ranks
        fD;             % set of dims for fibres
        bI = 1;         % base index
        fI = 2;         % fibre index
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fdi;             % fibre dictionary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoFB(varargin)
            %%%%%%%%%%%%%%%%
            % init calls as [baseRank,fibreDims]
            %   example - phytoFB(1,3) is a curve in three space
            %                   a curve with the fibre structure of a 
            %                   three dim vector space
            %   example - phytoFB(1,[3 3]) is a curve with a second rank
            %                   tensor attach to the curve. the tensor is
            %                   R^3 tensor R^2
            %   eample - phytoFB(2,3) is a surface in three space. the
            %                   along the surface is a vector fibre.
            obj.fdi = dictionary();
            %%%%%%%%%%%%%%%%
            % init dictionary
            if nargin == 2
                % set base rank, base fibre rank, base fibre dims
                obj.bR = varargin{1};
                obj.fR{1} = numel(varargin{2});
                obj.fD{1} = varargin{2};
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of base spaces
        function [n] = nBase(obj)            
            n = size(obj.g0,obj.bI);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of fibre spaces
        function [n] = nFibre(obj)            
            n = size(obj.g0,obj.fI);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % base rank  
        function [n] = baseRank(obj)            
            n = obj.bR;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % base dims
        function [n] = baseDims(obj)
            n = NaN;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get rank of fibre
        function [n] = fibreRank(obj,fI)
            if nargin == 1;fI = 1;end
            n = obj.fR{fI};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put trial
        function [] = putTrial(obj,T,bT,fT)
            %%%%%%%%%%%%%%%%%%%%%%%%
            % set default trial and fibre index
            if nargin == 2;bT = nBase(obj)+1;fT = 1;end
            % set default fibre index
            if nargin == 3;fT = 1;end
            %%%%%%%%%%%%%%%%%%%%%%%%
            % bundle split
            T = phytoFB.bfSplit(T,obj.bR);
            %%%%%%%%%%%%%%%%%%%%%%%%
            % store the data point
            for e = 1:size(T,1)
                obj.g0(bT,fT).b(e).f = T(e,:);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get Trial
        function [t] = getTrial(obj,bT,fT)
            % set default trial and fibre index
            if nargin == 2;fT = 1;end
            % hardcoded problem
            t = obj.g0(bT,fT);
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
        % used to split base-fibre complex into base-fibre stream.
        % used because all base and fibres are stored as stream form.
        function d = bfStream(d,bR)
            sz = size(d);
            szN = [prod(sz(1:bR)) prod(sz(bR+1:end))];
            d = reshape(d,szN);
        end
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % used to split base-fibre complex into base-fibre stream.
        % used because all base and fibres are stored as stream form.
        function d = bfStream(d,sz)                        
            d = reshape(d,szN);
        end
        %}
    end
    
    methods (Abstract)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % implement normalize
        %[] = normalize(obj,degree,trialI,fibreI);
    end
end