classdef myM < matlab.mixin.Copyable
    % my manifold class
    properties
        %%%%%%%%%%%%%%
        d;      % domain tensor
        c;      % codomain tensor
        %%%%%%%%%%%%%%
        f;      % protoypical set of functionals on tensor space
        mf;     % homogenous yet not parallel transported
                % distribution of functionals on the manifold
        cf;     % set of clustering functionals on the whole of the domain
        %%%%%%%%%%%%%%
        bv;     % set of tensor basis on manifold - for domain
        cbv;    % set of cotensors on manifold - for codomain
        %%%%%%%%%%%%%%
    end
    
    methods
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [obj] = myM(varargin)
                % init domain and codomain to []
                obj.d = [];
                obj.c = [];
                % init prototypical functional space
                obj.f = myF();
                obj.cf = myF();
                % init domain and codomain
                if nargin >= 1
                    setDomain(obj,varargin{1});
                end
                if nargin >= 2
                    setcoDomain(obj,varargin{2});
                end
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set domain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = setDomain(obj,d)
                if isa(d,'myT')
                    obj.d = d;
                else
                    obj.d = myT(d);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % clear domain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = clearDomain(obj)
                obj.d = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % clear co-domain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = clearcDomain(obj)
                obj.c = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % clear data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = clearData(obj)
                obj.clearDomain();
                obj.clearcDomain();
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set codomain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = setcoDomain(obj,c)
               if isa(c,'myT')
                    obj.c = c;
                else
                    obj.c = myT(c);
                end
            end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % add cluster method @ loc
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = addClusterMethod(obj,meth,loc)
                if nargin == 2
                    obj.cf.addFunctional(meth);
                end
                if nargin == 3
                    obj.cf.addFunctional(meth,loc);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call to cluster
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [obj] = cluster(obj)
                % vectorize the domain tensor
                obj.d.vec();
                % eval functional set on tensor
                obj.cf.i(obj.d,obj.c);
                % invert vec operation
                obj.d.i();
            end        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % map tensor set to centers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [ret] = mapTocenter(obj,d)                
                % vectorize the domain                
                d.vec();
                % map over all manifolds
                ret = obj.cf.e(d);
                % invert
                d.i();
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % decompose
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = decompose(obj)
                obj.decomposeDomain();
                obj.decomposecDomain();
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % decompose domain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = decomposeDomain(obj)
                % group
                grp = obj.mapTocenter(obj.d);

                % loop over manifolds
                for manifold = 1:numel(grp)
                    % obtain centers            
                    UQ = unique(grp{manifold});
                    % loop over clusters
                    for u = 1:numel(UQ)
                        % get the subtensor index
                        idx = find(grp{manifold} == UQ(u));
                        % create sub tensor
                        sT = obj.d.subset(idx);
                        % decompose and basis @ P
                        obj.bv{manifold,UQ(u)} = sT.decompose();
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % decompose codomain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = decomposecDomain(obj)
                % group
                grp = obj.mapTocenter(obj.d);

                % loop over manifolds
                for manifold = 1:numel(grp)
                    % obtain centers            
                    UQ = unique(grp{manifold});
                    % loop over clusters
                    for u = 1:numel(UQ)
                        % get the subtensor
                        idx = find(grp{manifold} == UQ(u));
                        % create sub tensor
                        sT = obj.c.subset(idx);
                        % decompose and basis @ P
                        obj.cbv{manifold,UQ(u)} = sT.decompose();
                    end
                end
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % add function on tensor space
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = addFunctional(obj,meth,loc)
                if nargin == 2
                    obj.f.addFunctional(meth);
                end
                if nargin == 3
                    obj.f.addFunctional(meth,loc);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % invert - fit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [] = i(obj)
                % group
                grp = obj.mapTocenter(obj.d);
                
                % loop over manifolds
                for manifold = 1:numel(grp)
                    % obtain centers            
                    UQ = unique(grp{manifold});
                    % loop over clusters
                    for u = 1:numel(UQ)
                        % current group number
                        grpNUM = UQ(u);
                        % distribute functionals
                        obj.mf{manifold,grpNUM} = copy(obj.f);
                        
                        % get the subtensor index
                        idx = find(grp{manifold} == grpNUM);
                        % create sub tensor for doamin
                        sTd_whole = obj.d.subset(idx);
                        sTc_whole = obj.c.subset(idx);
                        % project into basis vectors for both domain and                        
                        sTd_comp = obj.bv{manifold,grpNUM}.project(sTd_whole);
                        sTc_comp = obj.cbv{manifold,grpNUM}.project(sTc_whole);
                        
                        % create inverse mapping for subtensors
                        obj.mf{manifold,grpNUM}.i(sTd_comp,sTc_comp);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % eval
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [ret] = e(obj,d)                
                % group
                grp = obj.mapTocenter(d);
                ret = [];
                
                % loop over manifolds
                for manifold = 1:numel(grp)
                    % obtain centers            
                    UQ = unique(grp{manifold});
                    % loop over clusters to which the data vectors map to
                    for u = 1:numel(UQ)
                        % current group number
                        grpNUM = UQ(u);
                        % get the subtensor index of those which map to the
                        % uth center
                        idx = find(grp{manifold} == grpNUM);
                        % subset of data
                        sTd_whole = d.subset(idx);                        
                        % obtain the components
                        d_comp = obj.bv{manifold,grpNUM}.project(sTd_whole);                        
                        % eval the manifold function
                        cd_comp = obj.mf{manifold,grpNUM}.e(d_comp);                        
                        
                        
                        obj.cbv{manifold,grpNUM}.construct(cd_comp)
                        ret{manifold}(idx,:)
                    end
                end
            end
    end
end