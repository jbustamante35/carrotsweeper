classdef myT < matlab.mixin.Copyable & viewable & persistable
    % my tensor class
    properties
        d; % data
        h; % operation history        
    end
    % methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myT(varargin)
            if (nargin == 0)
                obj.d = [];
                obj.h = [];
            else
                if ~isjava(varargin{1})
                    obj.d = varargin{1};
                    obj.h = [];
                else
                    obj = myT.fromBson(varargin{1});
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % persist - NOT DONE
        function [] = persist(obj)
            % nothing yet
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % toBson
        function [xferO] = toBson(obj,xferO)
            import phytoG.locked.BdataObjects.geometry.implementations.*;    
            if (nargin == 1);xferO = phytoTensor();end            
            xferO.setData(mat2bson(obj.d));
            xferO.setProp('history',mat2bson(obj.h));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute
        function [] = permute(obj,p)
            % operate
            obj.d = permute(obj.d,p);
            % record history of permute
            obj.h(end+1).f = func2str(@(obj,x)ipermute(obj,x));
            obj.h(end).x = p;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ipermute
        function [] = ipermute(obj,p)
            % inverse permute
            obj.d = ipermute(obj.d,p);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape
        function [] = reshape(obj,nS)
            oldSZ = size(obj.d);
            obj.d = reshape(obj.d,nS);
            % record history
            obj.h(end+1).f = func2str(@(obj,x)ireshape(obj,x));
            obj.h(end).x = oldSZ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inverse reshape
        function [] = ireshape(obj,nS)
            obj.d = reshape(obj.d,nS);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute - remove later
        %{
        function [] = p(obj,x)
            permute(obj,x);
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inverse permute - removed
        %{
        function [] = ip(obj,v)
            obj.d = ipermute(obj.d,v);            
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = r(obj)
            % reshape
            x = size(obj.d);
            obj.d = reshape(obj.d,[x(1) prod(x(2:end))]);
            % record history
            obj.h(end+1).f = func2str(@(obj,x)ir(obj,x));
            obj.h(end).x = x;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = ir(obj,r)
            obj.d = reshape(obj.d,r);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute-reshape
        function [] = pr(obj,v)
            % v(i).d = ith grouping with dims
            % prod of ith grouping is size
            opand{1} = [];            
            opand{2} = size(obj.d);
            % stack the permute
            for i = 1:numel(v)
                opand{1} = [opand{1} v(i).d];
                sz(i) = prod(opand{2}(v(i).d));
            end
            % call to permute - with record history
            permute(obj,opand{1});
            % call to reshape - wit record history
            reshape(obj,sz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ipermute-reshape
        function [] = ipr(obj,v)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert dimension
        function [] = insertDim(obj,n)
            sz = size(obj.d);
            dims = ndims(obj.d);
            sz = [sz(1:(n-1)) 1 sz(n:dims)];
            obj.d = reshape(obj.d,sz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clear history
        function [] = clearHistory(obj)
            obj.h = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cycle permute
        function [] = cp(obj,k)
            % generate permute
            n = ndims(obj.d);
            p = 1:(n-1);
            p = circshift(p,[0 -k]);
            % perform permute
            obj.permute([p n]);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % roll back
        function [] = i(obj)
            % perform inverse(s)
            H = numel(obj.h);
            for e = 0:(H-1)
                index = H - e;
                func = str2func(obj.h(index).f);
                func(obj,obj.h(index).x);
            end
            % clean up
            obj.h = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rotate last rank to first (trials) to dim1
        function [] = vec(obj)
            % construct permute vector
            p = 1:ndims(obj.d);
            p = circshift(p,[0 1]);
            % permute
            obj.permute(p);
            % reshape
            obj.r();    
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate polynomial index
        function [] = gen_pi(obj)
            sz = size(obj.d);
            cumprod([1 sz]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % n-mode multiplication
        function [] = nm_mult(obj,m,n)
            % create dimension vector - switch order of [1,n]
            p = 1:ndims(obj.d);
            p(1) = n;
            p(n) = 1;
            % permute to focus on dim(n)
            obj.p(p);
            % reshape pre-mult
            obj.r();
            % perform multiplication
            obj.d = m*obj.d;
            % alter inverse history
            obj.h(end).x(1) = size(m,1);
            % inverse reshape
            obj.i();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % decompose
        function [bv] = decompose(obj)
            bv = ipca(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % subset via trial dim
        function [t] = subset(obj,idx)
            % vectorize
            obj.vec();
            %%%%%%%%%%%%%%%%
            % get subset
            t = obj.d(idx,:);
            % reshape subset tensor
            sz = obj.h(end).x;
            sz(1) = numel(idx);
            t = reshape(t,sz);
            t = ipermute(t,obj.h(1).x);
            % obtain tensor
            t = myT(t);
            %%%%%%%%%%%%%%%%
            % invert tensor
            obj.i();            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rank
        function [r] = rank(obj)
            r = ndims(obj.d);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number trials
        function [t] = nT(obj)
            t = size(obj.d,ndims(obj.d));
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view
        function [] = view(obj,h,vProps)
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);
            if ~isempty(obj.view_rep)
                tensorView(obj.view_rep,uProps,h);
            else
                tensorView(obj.d,uProps,h);
            end
        end
        % removed June 3, 2013
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view
        function [r] = rep(obj)
            r = obj.d;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsref
        function [r] = subsref(obj,S)
            if strcmp(S(1).type,'.')
                % put note when this happens: when nargout == 0 and error
                % re return value
                if nargout == 0
                    % ans: when there is a function call and S is numel==2
                    builtin('subsref',obj,S)
                else
                    r = builtin('subsref',obj,S);
                end
            else
                r = subsref(obj.d,S);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [obj] = subsasgn(obj,S,B)
            if strcmp(S(1).type,'.')                
                builtin('subsasgn',obj,S);                
            else
                obj.d = subsasgn(obj.d,S,B);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload plus
        function [r] = plus(obj,A)
            r = obj.copy();
            r.d = obj.d + A.d;
        end
    end
    % static methods
    methods (Static)        
        function [out] = fromBson(varargin)
            in = varargin{1};
            if (nargin == 1);out = myT();
            else out = varargin{2};end            
            out.d = bson2mat(in.getData());
            out.h = bson2mat(in.getProp('history'));
        end
        
    end
end

%{
% testing
data = rand(10,3,8,4);
T = myT(data);
V(1).d = [3 4];
V(2).d = [2 1];
T.pr(V);
T.i();
%}