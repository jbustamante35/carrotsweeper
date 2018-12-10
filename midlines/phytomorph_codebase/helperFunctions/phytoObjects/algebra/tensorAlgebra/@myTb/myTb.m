classdef myTb < myT & persistable
    % my tensor class
    properties
      bR;       % base rank
      fR;       % fibre rank
      tR;       % total rank
      sF = 1;   % shape flag
      oS;       % orginal size
      name;     % name of tensor bundle
    end
    
    % methods
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myTb(varargin)
            obj = obj@myT();
            if (nargin >= 2)
                % set the base rank, fibre rank, shapeFlag
                setAllRank(obj,varargin{1},varargin{2});           
                obj.sF = 1;
            end
            if (nargin == 3)
                setData(obj,varargin{3});
            end
            if (nargin == 1)
                if isjava(varargin{1})
                    obj = myTb.fromBson(varargin{1});
                end
            end
        end
        % toBson
        function [xferO] = toBson(obj,xferO)
            import phytoG.locked.BdataObjects.geometry.implementations.*;    
            if (nargin == 1);xferO = phytoTensor();end
            
            oT = toBson@myT(obj,xferO);
            oT.setProp('bR',mat2bson(obj.bR));
            oT.setProp('fR',mat2bson(obj.fR));
            oT.setProp('tR',mat2bson(obj.tR));
            oT.setProp('sF',mat2bson(obj.sF));
            oT.setProp('oS',mat2bson(obj.oS));
            oT.setProp('name',mat2bson(obj.name));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set base rank
        function [] = setAllRank(obj,bR,fR)
            setBaseRank(obj,bR);
            setFibreRank(obj,fR);
            obj.tR = bR + fR;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set base rank
        function [r] = setBaseRank(obj,r)
            obj.bR = r;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get base rank
        function [r] = getBaseRank(obj)
            r = obj.bR;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set fibre rank
        function [] = setFibreRank(obj,r)
            obj.fR = r;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get fibre rank
        function [r] = getFibreRank(obj)
            r = obj.fR;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get fibre rank
        function [r] = getTotalRank(obj)
            r = getBaseRank(obj) + getFibreRank(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % base size NEED TO UPDATE
        function [s] = baseSize(obj)
            s = 0;
            if obj.fR > 0
                % get orginal size of total space
                orgSize = totalSize(obj);
                % base Rank
                bR = getBaseRank(obj);
                % get fibreSize
                s = orgSize(1:bR);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fibre size
        function [s] = fibreSize(obj)
            s = 0;
            if obj.fR > 0
                % get orginal size of total space
                orgSize = totalSize(obj);
                % base Rank
                bR = getBaseRank(obj);
                % get fibreSize
                s = orgSize(bR+1:end);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fibre size
        function [s] = totalSize(obj)
            % get orginal size of total space
            s = obj.oS;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set fibre bundle
        function [] = setData(obj,d)
            % set data
            obj.d = d;
            % set orginal size - pad with ones
            s = size(obj.d);            
            obj.oS = [s ones(1,obj.tR-numel(s))];
            % fold
            fold(obj,obj.bR);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fold the tensor into base and fibre
        function [] = fold(obj,i)
            sz = totalSize(obj);
            if i > 0
                sz = [prod(sz(1:i)) prod(sz(i+1:end))];
            end
            reshape(obj,sz);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsref
        function [r] = subsref(obj,S)
            % if . type
            if strcmp(S(1).type,'.')
                if any(strcmp(methods(obj),S(1).subs))
                    try r = builtin('subsref',obj,S(1:2));catch ME;builtin('subsref',obj,S(1:2));end
                    S(1:2) = [];
                else
                    try r = builtin('subsref',obj,S(1));catch ME;builtin('subsref',obj,S(1));end
                    S(1) = [];
                end
            % if () type    
            elseif strcmp(S(1).type,'()')
                try r = indexT(obj,S(1));catch ME;indexT(obj,S(1));end
                S(1) = [];
            % if {} type
            elseif strcmp(S(1).type,'{}')
                try r = builtin('subsref',obj,S(1));catch ME;builtin('subsref',obj,S(1));end
                S(1) = [];
            end
            % recursive call(s)
            if numel(S) > 1
               % call to the next index level
               try r = builtin('subsref',r,S);catch ME;builtin('subsref',r,S);end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [obj] = subsasgn(obj,S,F)
            % generate index
            v = genI(obj,S);
            % reshape fibre
            F = myTb.reshapeFibre(F);
            % subs
            obj.d(v,:) = F;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % distribute fibre-wise function over base
        function [r] = distrib(obj,f)
            r = myHS_X('*');
            for e = 1:size(obj.d,1)
                ele = f(obj.d(e,:));
                r.putElement(ele);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % index Tensor
        function [v] = vectorize(obj,v)
            v = [v obj.d(:)'];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % helpers - private
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate polynomial index
        function [s] = gen_pi(obj)
            s = totalSize(obj);
            s = cumprod([s]);
            s = [1 s];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % index Tensor
        function [r] = indexT(obj,S)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the data from the block
            % if S has index values            
            if ~isempty(S.subs)
                
                if numel(S.subs) ~= 1
                    % generate index
                    v = genI(obj,S);
                else
                    v = S.subs{1};
                end
                
                % look up in data block
                r = obj.d(v,:);
            % if there is no index than
            % baseR == 0
            else
                r = obj.d;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reshape
            if obj.sF
                % if the base and fibre rank are not zero
                if getBaseRank(obj) ~= 0 & getFibreRank(obj) ~= 0
                    nsz = [size(r,1) fibreSize(obj)];
                    % if resize is needed - only if fibreSize > 1
                    if numel(nsz) > 1
                        r = reshape(r,nsz);
                    end                     
                % if the baserank is zero
                elseif getBaseRank(obj) == 0
                    nsz = fibreSize(obj);
                    if getFibreRank(obj) ~= 0
                        r = reshape(r,nsz);
                    end
                % if the fibre rank is zero
                elseif getFibreRank(obj) == 0
                    r = r;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate index
        function [v] = genI(obj,S)
            v = 0;
            if ~isempty(S.subs)
                % if : as index
                v = colon(obj,S(1).subs{end},numel(S(1).subs));
                % raster S to index
                for e = 1:numel(S(1).subs)-1
                    curV = colon(obj,S(1).subs{end-e},numel(S(1).subs)-e);
                    mag = ones(1,numel(curV));
                    v = kron(v,mag);
                    v = [repmat(curV,[1 size(v,2)/size(curV,2)]);v];
                end
                % generate poly
                s = gen_pi(obj);
                % get poly part
                s = s(1:size(v,1));
                % generate ind
                v = s*(v-1)+1;
            end 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload end
        function [r] = end(obj,k,T)
            r = obj.h(end).x(k);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload colon
        function [v] = colon(obj,v,e)
            % if colon as index then stack 1:N
            % where N is the dims of the vector
            % and N is stored in the history of the 
            % reshape
            if isa(v,'char')
                if strcmp(v,':')
                    v = 1:obj.h(end).x(e);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % helpers - private
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload plus
        function [r] = plus(obj,A)
            r = obj.copy();
            r.d = obj.d + A.d;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload plus
        function [r] = minus(obj,A)
            r = obj.copy();
            r.d = obj.d - A.d;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload mult
        function [r] = mtimes(obj,A,n)
            if nargin == 2
                r = A.copy();
                r.d = obj.d*A.d;
                r.view_rep = obj.d*A.view_rep;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape Fibre
        function [F] = reshapeFibre(F,sz)
            % reshape B
            sz = size(F);
            sz = [sz(1) prod(sz(2:end))];
            F = reshape(F,sz);
        end        
        % fromBson
        function [out] = fromBson(varargin)
            in = varargin{1};            
            
            if (nargin == 1);out = myTb();
            else;out = varargin{2};end
            
            out = myT.fromBson(in,out);
            out.bR = bson2mat(in.getProp('bR'));
            out.fR = bson2mat(in.getProp('fR'));
            out.tR = bson2mat(in.getProp('tR'));
            out.sF = bson2mat(in.getProp('sF'));
            out.oS = bson2mat(in.getProp('oS'));
            out.name = bson2mat(in.getProp('name'));
        end
    end
end