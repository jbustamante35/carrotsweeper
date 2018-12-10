classdef myT < matlab.mixin.Copyable & viewable
    % my tensor class
    properties
        d; % data
        h; % operation history
        %v; % view type and para
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
                obj.d = varargin{1};
                obj.h = [];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute
        function [] = p(obj,x)
            obj.d = permute(obj.d,x);
            % record history
            obj.h(end+1).f = @(x)ip(obj,x);            
            obj.h(end).x = x;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inverse permute
        function [] = ip(obj,v)
            obj.d = ipermute(obj.d,v);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape
        function [] = r(obj)
            x = size(obj.d);
            obj.d = reshape(obj.d,[x(1) prod(x(2:end))]);
            % record history
            obj.h(end+1).f = @(x)ir(obj,x);
            obj.h(end).x = x;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inverse reshape
        function [] = ir(obj,r)
            obj.d = reshape(obj.d,r);
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
        function [] = clearHistory()
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
            obj.p([p n]);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % roll back
        function [] = i(obj)
            % perform inverse(s)
            H = numel(obj.h);
            for e = 0:(H-1)
                index = H - e;
                obj.h(index).f(obj.h(index).x);
            end
            % clean up
            obj.h = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rotate trials to dim1
        function [t] = vec(obj)
            % construct permute vector
            p = 1:ndims(obj.d);
            p = circshift(p,[0 1]);
            
            % permute
            obj.p(p);
            
            % reshape
            obj.r();    
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape
        function [] = reshape(obj,nS)
            oldSZ = size(obj.d);
            obj.d = reshape(obj.d,nS);
            % record history
            obj.h(end+1).f = @(x)ireshape(obj,x);
            obj.h(end).x = oldSZ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inverse reshape
        function [] = ireshape(obj,nS)
            obj.d = reshape(obj.d,nS);
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
        function [] = view(obj,h,frame,vProps)
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);
            if ~isempty(obj.view_rep)
                tensorView(obj.view_rep,uProps,h);
            else
                tensorView(obj.d,uProps,h);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view
        function [r] = rep(obj,frame,type)
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
    end
end