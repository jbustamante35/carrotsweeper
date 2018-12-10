classdef stensor < tensor & hasSignature
    properties
        
    end
    
    methods
        % constructor
        function [obj] = tensor(M,sig)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting construction operation:\n']);tm = clock;
            if nargin == 0
                M = [];
            end
            % if the first object isa tensor -  then init properly
            if isa(M,'tensor')
                obj.M = M.M;
                obj.signature = M.signature;
            else
                % remove the singleton dims
                M = squeeze(M);
                % transpose/permute - rows and column vectors are not things in this universe
                if size(M,1) == 1
                    M = permute(M,[2 1]);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % calculate needed values
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nonTypedIndexNumber = sum(size(M) == 1);
                potentialTypedIndexNumber = ndims(M);
                typedIndexNumber = potentialTypedIndexNumber - nonTypedIndexNumber;
                if nargin == 1 | nargin == 0
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % assumed signature
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        for idx = 1:typedIndexNumber
                            obj.signature(idx) = TE(tE());
                        end
                else
                    if ~isa(sig,'TE')
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % assumed signature
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if numel(sig) == 1
                            sig = repmat(sig,[1 typedIndexNumber]);
                        end
                        for idx = 1:typedIndexNumber
                            obj.signature(idx) = TE(tE(sig(idx)));
                        end
                    else
                        obj.signature = sig; 
                    end

                end
                if ~isa(M,'double')
                    M = double(M);
                end
                obj.M = M;
            end
            % record the operation string
            opString = ['<-construct:'];
            postStr = obj.sig2str;
            fprintf([postStr '=' opString '.\n']);
            fprintf(['\tending construction operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
        
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'{}')
                 fprintf(['error in index-{} for tensor.\n']);
            elseif strcmp(S(1).type,'.')
                try
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                catch ME
                    getReport(ME);
                    fprintf(['error in index-. for tensor.\n']);
                end
            elseif strcmp(S(1).type,'()')
                [varargout{1:nargout}] = builtin('subsref',obj.M,S);
            end
        end
        
        % get the size of the tensor
        function [r] = size(obj)
            r = size(obj.M);
        end
        
        % get the order of the tensor
        function [r] = order(obj)
            r = size(obj.M);
            r(r==1) = [];
            r = numel(r);
        end
        
        % return the signature of the tensor as array of tE
        function [r] = sig(obj,idx)
            if nargin == 1
                idx = 1:numel(obj.signature);
            end
            r = [];
            strg = '';
            for e = idx
                r = [r obj.signature(e).Byte];
            end
        end
        
        % convert the signature to string
        function [r] = sig2str(obj)
            r = toString(obj.signature);
        end
        
        
        % permute the tensor
        function [] = permute(obj,order)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting permute operation:\n']);tm = clock;
            % record the operation string
            opString = ['<-p(' strrep(num2str(order),'  ',',') ')|'];
            % get the signature string pre operation
            preStr = obj.sig2str;
            % apply the permutation for the tensor
            if numel(order) == 1
                fixedOrder = [order 2];
            else
                fixedOrder = order;
            end
            obj.M = permute(obj.M,fixedOrder);
            % apply the permatation for the signature
            obj.signature = obj.signature.swap(order);
            % get the signature string post operation
            postStr = obj.sig2str;
            fprintf([postStr '=' opString preStr '|>.\n']);
            fprintf(['\tending permute operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
        
        % reshape the tensor object
        function [] = reshape(obj,newSZ)
           obj.M = reshape(obj.M,newSZ);
        end
        
        % fold the tensor object
        function [] = fold(obj,foldDirections)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting fold operation:\n']);tm = clock;
            % construct the string for operation
            r = '{';
            for l = 1:numel(foldDirections)
                r = [r '[' strrep(num2str(foldDirections{l}),'  ',',') '],'];
            end
            r(end) = [];
            r = [r '}'];
            opString = ['<-f(' r ')|'];
            % get the preoperation signature
            preStr = obj.sig2str();
            
            
            
            % decode the fold instructions into a permute instruction operand
            pVec = [];
            oldSZ = obj.size();
            newSZ = zeros(1,numel(foldDirections));
            for e = 1:numel(foldDirections)
                pVec = [pVec foldDirections{e}];
                newSZ(e) = [prod(oldSZ(foldDirections{e}))];
            end
            % if the newSZ is a column/row vector then make [column 1]
            if numel(newSZ) == 1
                newSZ = [newSZ 1];
            end
            
            
            % construct the reshape instruction operad for post permute
            str = 1;
            for e = 1:numel(foldDirections)
                stp = (str+numel(foldDirections{e})-1);
                sigDirections{e} = str:stp;
                str = stp + 1;
            end
            
            % check that all dims are accounted for
            dF = setdiff(min(pVec):max(pVec),1:numel(oldSZ));
            if isempty(dF)
                obj.permute(pVec);
                
                obj.reshape(newSZ);
            else
                fprintf(['error in fold operation due to lack of dims accounting.\n']);
            end
            
            % fold the signature
            obj.signature = obj.signature.fold(sigDirections);
            
            % get the signature post operation
            postStr = obj.sig2str();
            fprintf([postStr '=' opString '||' preStr '|>.\n']);
            fprintf(['\tending fold operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
        
        function [r] = sum(obj,order)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting sum operation:\n']);tm = clock;
            opString = ['||(' '*' ')x(' num2str(order) ')||'];
            preStr = obj.sig2str;
            r = tensor(sum(obj.M,order),obj.signature.removeOrder(order));
            postStr = r.sig2str;
            fprintf([postStr '<-' opString preStr '|>.\n']);
            fprintf(['\tending sum operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
        
        function [r] = prod(obj,order)
            r = tensor(prod(obj.M,order),obj.signature.removeOrder(order));
        end
        
        function [] = applyAlong(obj,func,order)
            obj.M = func(obj.M,order);
            obj.signature.removeOrder(order);
        end
        
        function [r] = or(objA,objB)
            if objB.order == 2 && objA.order == 1
                r = tensor.dotProduct(objA,objB,[1 1]);
            elseif objA.order == 2 && objB.order == 1
                r = tensor.dotProduct(objA,objB,[2 1]);
            else
                r = tensor.dotProduct(objA,objB,[1 1]);
            end
           
        end
        
        function [r] = gt(objA,objB)
            objA.fold({[1:objA.order()]});
            r = objA;
        end
        
        function [r] = lt(objA,objB)
            objB.fold({[1:objB.order()]});
            r = objB;
        end
        
        
        
        function [r] = plus(objA,objB)
            if isa(objA,'double')
                r = tensor(objB.M + objA,objB.signature);
            elseif isa(objB,'double')
                r = tensor(objA.M + objB,objA.signature);
            else
                if tensor.orderCheck(objA,objB)
                    if tensor.dimsCheck(objA,objB)
                        for e = 1:numel(objA.order)
                            c(e) = tensor.typeCheck(objA,objB,[e e]);
                        end
                        if all(c)
                            r = tensor(objA.M+objB.M,objA.signature);
                        end
                    end
                end
            end
                
                
        end
        
        function [r] = getExpectedVector(obj,d)
            sz = size(obj.M);
            v = ones(sz(d),1)/sz(d);
            [~,sig,~,~] = pullSignature(obj,d);
            r = tensor(v,sparse(~sig(:,1)));
        end
    end
    
    
    methods (Static)
        
        function [c] = typeCheck(a,b,d)
            c = a.sig(d(1)) == b.sig(d(2));
        end
        
        function [c] = dimsCheck(a,b,d)
            if nargin == 2
                if tensorContainer.orderCheck(a,b)
                    d = [1:a.order() ; 1:b.order()]';
                else
                    fprintf(['order does not match therefore dims can not!\n']);
                    c = false;
                end
            end
            szA = size(a);
            szB = size(b);
            
           
            ca = szA(d(:,1));
            cb = szB(d(:,2));
            c = all(all(ca == cb,2),1);
        end
        
        function [c] = orderCheck(a,b)
            c = a.order() == b.order();
        end
        
        % dirty due to lack of checks on dims,types,order
        function [r] = bsxfun(binaryOp,objA,objB,shiftInstructions)
            a = shiftdim(objA.M,shiftInstructions(1));
            b = shiftdim(objB.M,shiftInstructions(2));
            r = bsxfun(binaryOp,a,b);
            r = tensor(r,objA.signature);
        end
        
        function [r] = dotProduct(a,b,d)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting inner-product operation:\n']);tm = clock;
            opString = ['||(' num2str(d(1)) ')x(' num2str(d(2)) ')||'];
            if isa(a,'tensor') && isa(b,'tensor')
                if tensor.dimsCheck(a,b,d)
                    [typeC] = tensor.typeCheck(a,b,d);
                    if typeC
                        astr = a.sig2str();
                        bstr = b.sig2str();
                        at = a.signature;
                        bt = b.signature;
                        a = a.M;
                        b = b.M;
                        r = TA.dot(a,b,d);
                        at = at.removeOrder(d(1));
                        bt = bt.removeOrder(d(2));
                        newt = [at bt];
                        if isempty(newt)
                            newt = TE(tE(0));
                        end
                        r = tensor(r,newt);
                    end

                end
                fprintf([r.sig2str() '= <|' astr opString bstr '|>.\n']);
            end
            fprintf(['\tending inner-product operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
        
        function [r] = tensorProduct(a,b)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting tensor-product operation:\n']);tm = clock;
            opString = ['||*||'];
            astr = a.sig2str();
            bstr = b.sig2str();
            r = TA.tensor(a.M,b.M);
            rSignature = [a.signature b.signature];
            r = tensor(r,rSignature);
            fprintf([r.sig2str() '= <|' astr opString bstr '|>.\n']);
            fprintf(['\tending tensor-product operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
    end
end
%{

%}






