classdef myE < matlab.mixin.Copyable & viewable
    % my basis vector class
    properties
        d;      % set of basis functions for each vector space
        bv;     % expanded form
        info;   % information about projection        
        default_threshold = .90;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myE(varargin)
            obj.d = [];
            obj.bv = [];
            obj.info = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add basis vector set
        function [obj] = addBasis(obj,e,u,n)
            if nargin <= 3
                n = numel(obj.d)+1;
            end
            obj.d(n).e = e;
            obj.d(n).u = u;
            obj.d(n).nv = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate basis functions
        function [] = generateBasisTensors(obj)
            % number of components
            nvec = fliplr(obj.rank());
            %%%%%%%%%%%%%%%%%%%%%
            % DEFAULTS 
            %%%%%%%%%%%%%%%%%%%%%
            if isempty(nvec)
                obj.selectBasisNumber();                 
                nvec = fliplr(obj.rank());
            end
            %%%%%%%%%%%%%%%%%%%%%
            % DEFAULTS
            %%%%%%%%%%%%%%%%%%%%%
            % prod nvec
            mNum = prod(nvec);    
            % loop over all combinations of basis vectors 
            % from each space
            for jmp = 1:mNum
                % get the coeffs for the forLoop
                ct = nCm((jmp-1),nvec); 
                % clear tensor basis
                tmp_bt = [];
                % generate basis tensors
                for space = 1:numel(ct)
                    v = obj.d(space).e(:,ct(space)+1);
                    
                    
                    
                    if space == 1
                        tmp_bt.s = numel(v);
                        tmp_bt.d = v;
                    else
                        tmp_bt = multv(tmp_bt,v);
                    end
                end        
                % reshape basis tensor                
                obj.bv = [obj.bv tmp_bt.d];
            end
            obj.bv = hilbertNormalize(obj.bv')';
            
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get R-rank
        function [v] = rank(obj)
            v = [];            
            for e = 1:numel(obj.d)
                if ~isempty(obj.d(e).nv)
                    v(e) = obj.d(e).nv;
                else
                    return;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get span
        function [v] = span(obj)
            v = [];            
            for e = 1:numel(obj.d)                
                v(e) = size(obj.d(e).e,1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get information
        function [] = attachInformation(obj,info)
           obj.info = info;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % selectBasisNumber
        function [] = selectBasisNumber(obj,threshold)            
            if nargin == 1
                threshold = obj.default_threshold;
            end
            
            for e = 1:numel(obj.d)
                fidx = find(obj.info(e).pe < threshold);
                if isempty(fidx)
                    fidx = 1;
                end
                obj.d(e).nv = fidx(end);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % project into
        function [comp] = project(obj,data)
            % generate basis tensors if needed
            if isempty(obj.bv)
                obj.generateBasisTensors();
            end
            % rotate the trials to front
            data.vec();
            % project components
            comp = data.d*obj.bv;
            % inverse
            data.i();
            % return components in tensor structure
            comp = myT(comp);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % project out-of
        function [sim] = construct(obj,comp)
            sim = obj.bv*comp.d';
            sim = reshape(sim,[obj.span size(sim,2)]);
            sim = myT(sim);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view
        function [] = view(obj,h,varargin)
        
        end
    end
end