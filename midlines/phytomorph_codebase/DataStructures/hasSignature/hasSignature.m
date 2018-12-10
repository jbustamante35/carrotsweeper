classdef hasSignature < handle
    properties
        signature = TE(tE());
    end
    
    methods
        % constructor
        function [obj] = hasSignature(signatureOrder,initSig)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting signature construction operation:\n']);tm = clock;
            if isa(initSig,'hasSignature')
                obj = initSig;
            elseif nargin == 1
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % assumed signature
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for idx = 1:signatureOrder
                        obj.signature(idx) = TE(tE());
                    end
            elseif nargin == 2
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % assumed signature
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if numel(initSig) == 1
                        initSig = repmat(initSig,[1 signatureOrder]);
                    end
                    for idx = 1:signatureOrder
                        obj.signature(idx) = TE(tE(initSig(idx)));
                    end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % record the operation string
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            opString = ['<-construct:'];
            postStr = obj.sig2str;
            fprintf([postStr '=' opString '.\n']);
            fprintf(['\tending  signature construction operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
            
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
        
        % permute the signature
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
           % not yet needed to signature - folding only
        end
        
        % fold the tensor object
        function [] = fold(obj,foldDirections)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting signature fold operation:\n']);tm = clock;
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
            
            
            % fold the signature
            obj.signature = obj.signature.fold(foldDirections);
            
            % get the signature post operation
            postStr = obj.sig2str();
            fprintf([postStr '=' opString '||' preStr '|>.\n']);
            fprintf(['\tending fold operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
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






